#!/usr/bin/env python3

import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, TextIO

from tqdm import tqdm


def open_text(path: str, mode: str = "rt") -> TextIO:
    """Open plain text or gzipped text based on filename suffix."""
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def load_rsid_set(path: str) -> Set[str]:
    """Load rsIDs (first token per line) into a set."""
    with open_text(path, "rt") as f:
        return {line.strip().split()[0] for line in f if line.strip()}


def safe_float(x: str) -> Optional[float]:
    """Parse float safely; return None if missing or invalid."""
    if x in (".", "", "NA", "nan", "NaN", "None", "NULL"):
        return None
    try:
        return float(x)
    except Exception:
        return None


def parse_format_indices(format_str: str, wanted: List[str]) -> Dict[str, int]:
    """
    Return mapping from FORMAT field name to index.
    Raise ValueError if any required field is missing.
    """
    keys = format_str.split(":")
    idx = {k: i for i, k in enumerate(keys)}
    missing = [k for k in wanted if k not in idx]
    if missing:
        raise ValueError(
            "VCF FORMAT missing fields: {}. Found: {}".format(missing, keys)
        )
    return idx


def count_non_header_lines(path: str) -> int:
    """Count non-header VCF lines for a real tqdm progress bar."""
    n = 0
    with open_text(path, "rt") as f:
        for line in f:
            if line and not line.startswith("#"):
                n += 1
    return n


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert IEU/OpenGWAS VCF(.vcf/.vcf.gz) to LDSC-compatible TSV(.tsv/.tsv.gz). "
            "Extracts ES, SE, LP, optional SS, and can optionally filter by rsID list."
        )
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input IEU/OpenGWAS VCF (.vcf or .vcf.gz)"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output LDSC TSV (.tsv or .tsv.gz)"
    )
    parser.add_argument(
        "--rsid-list",
        default=None,
        help="Optional rsID list (one per line) used for filtering"
    )
    parser.add_argument(
        "--rsid-col-name",
        default="SNP",
        help="Header name for rsID column in output (default: SNP)"
    )
    parser.add_argument(
        "--keep-only-rsids",
        action="store_true",
        help="Keep only variants present in --rsid-list"
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bar"
    )

    parser.add_argument(
        "--N",
        type=float,
        default=None,
        help="Constant sample size used when SS is missing in VCF FORMAT"
    )

    parser.add_argument(
        "--field-es",
        default="ES",
        help="FORMAT field name for effect size (default: ES)"
    )
    parser.add_argument(
        "--field-se",
        default="SE",
        help="FORMAT field name for standard error (default: SE)"
    )
    parser.add_argument(
        "--field-lp",
        default="LP",
        help="FORMAT field name for -log10(P) (default: LP)"
    )
    parser.add_argument(
        "--field-ss",
        default="SS",
        help="FORMAT field name for sample size (default: SS)"
    )

    args = parser.parse_args()

    vcf_path = args.vcf
    out_path = args.out

    if not Path(vcf_path).exists():
        raise FileNotFoundError("Input VCF not found: {}".format(vcf_path))

    out_parent = Path(out_path).parent
    out_parent.mkdir(parents=True, exist_ok=True)

    if args.N is not None and args.N <= 0:
        raise ValueError("--N must be > 0")

    rsids: Optional[Set[str]] = None
    if args.rsid_list:
        rsids = load_rsid_set(args.rsid_list)
        print(
            "Loaded {:,} rsIDs from {}".format(len(rsids), args.rsid_list),
            file=sys.stderr
        )
        if args.keep_only_rsids:
            print(
                "Filtering enabled: keeping only variants present in rsID list.",
                file=sys.stderr
            )

    wanted_base = [args.field_es, args.field_se, args.field_lp]
    wanted_with_ss = wanted_base + [args.field_ss]

    total = 0
    kept = 0
    removed = 0
    skipped_badline = 0
    skipped_missing = 0
    used_ss_count = 0
    used_constN_count = 0

    fmt_cache: Dict[str, Dict[str, int]] = {}

    progress_total = None
    if not args.no_progress:
        print("Counting variants for progress bar...", file=sys.stderr)
        progress_total = count_non_header_lines(vcf_path)
        print(
            "Found {:,} non-header variants in input VCF.".format(progress_total),
            file=sys.stderr
        )

    with open_text(vcf_path, "rt") as fin, open_text(out_path, "wt") as fout:
        fout.write(
            "\t".join([args.rsid_col_name, "A1", "A2", "BETA", "SE", "P", "N"]) + "\n"
        )

        iterator = fin
        if not args.no_progress:
            iterator = tqdm(
                fin,
                total=progress_total,
                desc="Reading VCF",
                unit="lines",
                dynamic_ncols=True
            )

        for line in iterator:
            if not line or line.startswith("#"):
                continue

            total += 1
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 10:
                skipped_badline += 1
                continue

            chrom = parts[0]
            pos = parts[1]
            vid = parts[2]
            ref = parts[3]
            alt = parts[4]
            fmt = parts[8]
            sample = parts[9]

            if vid in (".", ""):
                skipped_missing += 1
                continue

            if "," in alt:
                skipped_missing += 1
                continue

            if rsids is not None and args.keep_only_rsids and vid not in rsids:
                removed += 1
                continue

            if fmt not in fmt_cache:
                keys = fmt.split(":")
                has_ss = args.field_ss in keys
                wanted = wanted_with_ss if has_ss else wanted_base
                idx = parse_format_indices(fmt, wanted)
                idx["__has_ss__"] = 1 if has_ss else 0
                fmt_cache[fmt] = idx

            fmt_idx = fmt_cache[fmt]
            vals = sample.split(":")

            es = safe_float(vals[fmt_idx[args.field_es]]) if fmt_idx[args.field_es] < len(vals) else None
            se = safe_float(vals[fmt_idx[args.field_se]]) if fmt_idx[args.field_se] < len(vals) else None
            lp = safe_float(vals[fmt_idx[args.field_lp]]) if fmt_idx[args.field_lp] < len(vals) else None

            n_val: Optional[float] = None
            if fmt_idx.get("__has_ss__", 0) == 1:
                ss = safe_float(vals[fmt_idx[args.field_ss]]) if fmt_idx[args.field_ss] < len(vals) else None
                n_val = ss
                if n_val is not None:
                    used_ss_count += 1
            else:
                n_val = args.N
                if n_val is not None:
                    used_constN_count += 1

            if es is None or se is None or lp is None or n_val is None:
                skipped_missing += 1
                continue

            if n_val <= 0:
                skipped_missing += 1
                continue

            try:
                p = 10 ** (-lp)
            except OverflowError:
                skipped_missing += 1
                continue

            # OpenGWAS VCF convention:
            # ALT is treated as effect allele, REF as non-effect allele.
            a1 = alt
            a2 = ref

            fout.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    vid,
                    a1,
                    a2,
                    es,
                    se,
                    p,
                    int(round(n_val))
                )
            )
            kept += 1

    if kept == 0 and used_ss_count == 0 and args.N is None:
        print(
            "\nERROR: No variants were kept and SS was not present in FORMAT. "
            "Provide --N <sample_size>.\n",
            file=sys.stderr
        )

    print("\nDone.", file=sys.stderr)
    print("Input VCF: {}".format(vcf_path), file=sys.stderr)
    print("Output TSV: {}".format(out_path), file=sys.stderr)
    print("Total variants read (non-header): {:,}".format(total), file=sys.stderr)

    if total > 0:
        pct_kept = 100.0 * kept / total
        pct_removed = 100.0 * removed / total
        pct_badline = 100.0 * skipped_badline / total
        pct_missing = 100.0 * skipped_missing / total
    else:
        pct_kept = pct_removed = pct_badline = pct_missing = 0.0

    print("Kept: {:,} ({:.2f}%)".format(kept, pct_kept), file=sys.stderr)

    if rsids is not None and args.keep_only_rsids:
        print(
            "Removed (not in rsID list): {:,} ({:.2f}%)".format(removed, pct_removed),
            file=sys.stderr
        )

    print(
        "Skipped (bad lines): {:,} ({:.2f}%)".format(skipped_badline, pct_badline),
        file=sys.stderr
    )
    print(
        "Skipped (missing/invalid ES/SE/LP/N): {:,} ({:.2f}%)".format(
            skipped_missing, pct_missing
        ),
        file=sys.stderr
    )
    print("Used N from SS field: {:,}".format(used_ss_count), file=sys.stderr)
    print("Used constant --N: {:,}".format(used_constN_count), file=sys.stderr)


if __name__ == "__main__":
    main()
