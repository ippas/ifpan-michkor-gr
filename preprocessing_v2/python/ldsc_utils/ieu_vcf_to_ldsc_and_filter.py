#!/usr/bin/env python3
import argparse
import gzip
import sys
from typing import Dict, List, Optional, Set, Tuple
from tqdm import tqdm


def open_text(path: str, mode: str = "rt"):
    """Open plain text or gzipped text based on filename suffix."""
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def load_rsid_set(path: str) -> Set[str]:
    """Load rsIDs (first token per line) into a set."""
    with open_text(path, "rt") as f:
        return set(line.strip().split()[0] for line in f if line.strip())


def safe_float(x: str) -> Optional[float]:
    """Parse float safely; return None if missing/invalid."""
    if x in (".", "", "NA", "nan", "NaN", "None", "NULL"):
        return None
    try:
        return float(x)
    except Exception:
        return None


def parse_format_indices(format_str: str, wanted: List[str]) -> Dict[str, int]:
    """
    Return dict mapping FORMAT keys -> index.
    Raises if any wanted field is missing.
    """
    keys = format_str.split(":")
    idx = {k: i for i, k in enumerate(keys)}
    missing = [k for k in wanted if k not in idx]
    if missing:
        raise ValueError(f"VCF FORMAT missing fields: {missing}. Found: {keys}")
    return idx


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Convert IEU/OpenGWAS VCF(.gz) to LDSC-compatible TSV(.gz) and optionally filter by rsID list.\n"
            "Supports FORMAT with ES/SE/LP and optional SS (sample size). If SS is absent, use --N constant."
        )
    )
    ap.add_argument("--vcf", required=True, help="Input IEU/OpenGWAS VCF(.vcf or .vcf.gz)")
    ap.add_argument("--out", required=True, help="Output LDSC TSV(.tsv or .tsv.gz)")
    ap.add_argument("--rsid-list", default=None, help="Optional rsID list (one per line) for filtering")
    ap.add_argument("--rsid-col-name", default="SNP", help="Header name for rsID column (default: SNP)")
    ap.add_argument("--keep-only-rsids", action="store_true",
                    help="If set, keep only variants present in --rsid-list")
    ap.add_argument("--no-progress", action="store_true", help="Disable progress bar")

    # Fallback constant sample size (N)
    ap.add_argument("--N", type=float, default=None,
                    help="Constant sample size to write into output column N if SS is missing in VCF FORMAT.")

    # Allow overriding FORMAT field names
    ap.add_argument("--field-es", default="ES", help="FORMAT field name for effect size (default: ES)")
    ap.add_argument("--field-se", default="SE", help="FORMAT field name for standard error (default: SE)")
    ap.add_argument("--field-lp", default="LP", help="FORMAT field name for -log10(P) (default: LP)")
    ap.add_argument("--field-ss", default="SS", help="FORMAT field name for sample size N (default: SS)")

    args = ap.parse_args()

    if args.N is not None and args.N <= 0:
        raise ValueError("--N must be > 0")

    # Load rsID filter set if provided
    rsids: Optional[Set[str]] = None
    if args.rsid_list:
        rsids = load_rsid_set(args.rsid_list)
        print(f"Loaded {len(rsids):,} rsIDs from {args.rsid_list}", file=sys.stderr)
        if args.keep_only_rsids:
            print("Filtering enabled: keeping only variants present in rsID list.", file=sys.stderr)

    # Fields required always
    wanted_base = [args.field_es, args.field_se, args.field_lp]
    wanted_with_ss = wanted_base + [args.field_ss]

    # Stats
    kept = 0
    removed = 0
    total = 0
    skipped_badline = 0
    skipped_missing = 0
    used_ss_count = 0
    used_constN_count = 0

    # Cache FORMAT parsing
    fmt_idx: Optional[Dict[str, int]] = None

    with open_text(args.vcf, "rt") as fin, open_text(args.out, "wt") as fout:
        # Write LDSC header
        fout.write("\t".join([args.rsid_col_name, "A1", "A2", "BETA", "SE", "P", "N"]) + "\n")

        it = fin
        if not args.no_progress:
            it = tqdm(fin, desc="Reading VCF", unit="lines")

        for line in it:
            if not line or line[0] == "#":
                continue

            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                skipped_badline += 1
                continue

            # VCF columns
            vid = parts[2]      # rsID typically
            ref = parts[3]
            alt = parts[4]
            fmt = parts[8]
            sample = parts[9]

            # rsID filtering
            if rsids is not None and args.keep_only_rsids and (vid not in rsids):
                removed += 1
                continue

            # Re-parse FORMAT indices if changed
            if fmt_idx is None or fmt != fmt_idx.get("__fmt__", None):
                keys = fmt.split(":")
                has_ss = args.field_ss in keys
                wanted = wanted_with_ss if has_ss else wanted_base
                idx = parse_format_indices(fmt, wanted)
                idx["__fmt__"] = fmt
                idx["__has_ss__"] = 1 if has_ss else 0
                fmt_idx = idx

            vals = sample.split(":")

            # Extract ES/SE/LP
            es = safe_float(vals[fmt_idx[args.field_es]]) if fmt_idx[args.field_es] < len(vals) else None
            se = safe_float(vals[fmt_idx[args.field_se]]) if fmt_idx[args.field_se] < len(vals) else None
            lp = safe_float(vals[fmt_idx[args.field_lp]]) if fmt_idx[args.field_lp] < len(vals) else None

            # Determine N (sample size)
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

            # Validate required values
            if es is None or se is None or lp is None or n_val is None:
                skipped_missing += 1
                continue

            # Compute P from LP = -log10(P)
            # P = 10^(-LP)
            p = 10 ** (-lp)

            # LDSC expects A1/A2; common convention: A1=ALT(effect allele), A2=REF
            fout.write(f"{vid}\t{alt}\t{ref}\t{es}\t{se}\t{p}\t{int(n_val)}\n")
            kept += 1

    # If SS never present and N not provided, everything will be skipped; help user
    if kept == 0 and used_ss_count == 0 and args.N is None:
        print(
            "\nERROR: No variants were kept and SS was not present in FORMAT. "
            "Provide --N <sample_size> to write constant N.\n",
            file=sys.stderr
        )

    print("\nDone.", file=sys.stderr)
    print(f"Total variants read (non-header): {total:,}", file=sys.stderr)

    if total > 0:
        pct_kept = 100 * kept / total
        pct_removed = 100 * removed / total
        pct_skipped_missing = 100 * skipped_missing / total
        pct_badline = 100 * skipped_badline / total
    else:
        pct_kept = pct_removed = pct_skipped_missing = pct_badline = 0.0

    print(f"Kept: {kept:,} ({pct_kept:.2f}%)", file=sys.stderr)

    if rsids is not None and args.keep_only_rsids:
        print(f"Removed (not in rsID list): {removed:,} ({pct_removed:.2f}%)", file=sys.stderr)

    print(f"Skipped (bad lines): {skipped_badline:,} ({pct_badline:.2f}%)", file=sys.stderr)
    print(f"Skipped (missing/invalid ES/SE/LP/N): {skipped_missing:,} ({pct_skipped_missing:.2f}%)", file=sys.stderr)

    print(f"Used N from SS field: {used_ss_count:,}", file=sys.stderr)
    print(f"Used constant --N: {used_constN_count:,}", file=sys.stderr)


if __name__ == "__main__":
    main()