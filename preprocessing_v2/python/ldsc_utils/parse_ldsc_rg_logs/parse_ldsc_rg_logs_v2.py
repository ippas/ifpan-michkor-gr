#!/usr/bin/env python3

import re
import csv
import argparse
from pathlib import Path
from typing import Dict, Optional, List
from tqdm import tqdm


############################
# regex helper definitions #
############################

FLOAT_RE = r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?"
INT_RE = r"\d+"

RE_BEGINNING_ANALYSIS = re.compile(r"^Beginning analysis at\s+(.*)$", re.MULTILINE)
RE_FINISHED_ANALYSIS = re.compile(r"^Analysis finished at\s+(.*)$", re.MULTILINE)
RE_TOTAL_TIME = re.compile(r"^Total time elapsed:\s*(" + FLOAT_RE + r")s$", re.MULTILINE)

RE_OUT = re.compile(r"^--out\s+(.+)$", re.MULTILINE)
RE_RG = re.compile(r"^--rg\s+(.+)$", re.MULTILINE)
RE_REF_LD = re.compile(r"^--ref-ld-chr\s+(.+)$", re.MULTILINE)
RE_W_LD = re.compile(r"^--w-ld-chr\s+(.+)$", re.MULTILINE)

RE_READ_SUMSTATS = re.compile(
    r"Reading summary statistics from\s+(.+?)\s+\.\.\.", re.MULTILINE
)

RE_READ_SUMSTATS_SNPS = re.compile(
    r"Read summary statistics for\s+(" + INT_RE + r") SNPs\.", re.MULTILINE
)

RE_REF_SNPS = re.compile(
    r"Read reference panel LD Scores for\s+(" + INT_RE + r") SNPs\.", re.MULTILINE
)

RE_WEIGHT_SNPS = re.compile(
    r"Read regression weight LD Scores for\s+(" + INT_RE + r") SNPs\.", re.MULTILINE
)

RE_AFTER_MERGE_REF = re.compile(
    r"After merging with reference panel LD,\s+(" + INT_RE + r") SNPs remain\.", re.MULTILINE
)

RE_AFTER_MERGE_REG = re.compile(
    r"After merging with regression SNP LD,\s+(" + INT_RE + r") SNPs remain\.", re.MULTILINE
)

RE_AFTER_MERGE_SUMSTATS = re.compile(
    r"After merging with summary statistics,\s+(" + INT_RE + r") SNPs remain\.", re.MULTILINE
)

RE_VALID_ALLELES = re.compile(
    r"(" + INT_RE + r") SNPs with valid alleles\.", re.MULTILINE
)

RE_COMPUTING_RG = re.compile(
    r"Computing rg for phenotype\s+(\d+)/(\d+)", re.MULTILINE
)

RE_PHENO1_H2 = re.compile(
    r"Heritability of phenotype 1\s*-+\s*"
    r"Total Observed scale h2:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)\s*"
    r"Lambda GC:\s*(" + FLOAT_RE + r")\s*"
    r"Mean Chi\^2:\s*(" + FLOAT_RE + r")\s*"
    r"Intercept:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)\s*"
    r"Ratio:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)",
    re.MULTILINE | re.DOTALL
)

RE_PHENO2_H2 = re.compile(
    r"Heritability of phenotype 2(?:/\d+)?\s*-+\s*"
    r"Total Observed scale h2:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)\s*"
    r"Lambda GC:\s*(" + FLOAT_RE + r")\s*"
    r"Mean Chi\^2:\s*(" + FLOAT_RE + r")\s*"
    r"Intercept:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)\s*"
    r"Ratio:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)",
    re.MULTILINE | re.DOTALL
)

RE_GENCOV = re.compile(
    r"Genetic Covariance\s*-+\s*"
    r"Total Observed scale gencov:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)\s*"
    r"Mean z1\*z2:\s*(" + FLOAT_RE + r")\s*"
    r"Intercept:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)",
    re.MULTILINE | re.DOTALL
)

RE_GENCOR = re.compile(
    r"Genetic Correlation\s*-+\s*"
    r"Genetic Correlation:\s*(" + FLOAT_RE + r")\s*\((" + FLOAT_RE + r")\)\s*"
    r"Z-score:\s*(" + FLOAT_RE + r")\s*"
    r"P:\s*(" + FLOAT_RE + r")",
    re.MULTILINE | re.DOTALL
)

RE_SUMMARY_ROW = re.compile(
    r"Summary of Genetic Correlation Results.*?\n"
    r".*?\n"
    r"\s*(\S+)\s+(\S+)\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")\s+(" + FLOAT_RE + r")",
    re.MULTILINE | re.DOTALL
)

RE_ERROR = re.compile(
    r"(Traceback \(most recent call last\):.*|"
    r"ValueError:.*|"
    r"Exception:.*|"
    r"ERROR:.*|"
    r"Error:.*|"
    r"FloatingPointError:.*|"
    r"ZeroDivisionError:.*|"
    r"LinAlgError:.*)",
    re.MULTILINE | re.DOTALL
)


###################
# helper functions #
###################

def safe_float(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def safe_int(x: Optional[str]) -> Optional[int]:
    if x is None:
        return None
    try:
        return int(x)
    except ValueError:
        return None


def basename_without_sumstats(path_str: Optional[str]) -> Optional[str]:
    if not path_str:
        return None

    name = Path(path_str).name

    suffixes = [
        ".sumstats.gz",
        ".sumstats",
        ".gz",
        ".tsv",
        ".txt",
    ]

    for suffix in suffixes:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break

    return name


def find_all(pattern: re.Pattern, text: str) -> List[str]:
    return pattern.findall(text)


def find_first(pattern: re.Pattern, text: str) -> Optional[str]:
    m = pattern.search(text)
    if not m:
        return None
    if len(m.groups()) == 0:
        return m.group(0)
    return m.group(1)


########################
# log parsing function #
########################

def parse_ldsc_log(log_path: Path) -> Dict[str, Optional[object]]:
    text = log_path.read_text(errors="replace")

    row: Dict[str, Optional[object]] = {
        "log_file": str(log_path),
        "log_filename": log_path.name,

        "status": "OK",
        "error_message": None,

        "analysis_start": None,
        "analysis_end": None,
        "total_time_sec": None,

        "out_prefix": None,
        "ref_ld_chr": None,
        "w_ld_chr": None,
        "rg_arg_raw": None,

        "p1_path": None,
        "p2_path": None,
        "p1_id": None,
        "p2_id": None,

        "summary_stats_files_n": None,
        "computing_rg_current": None,
        "computing_rg_total": None,

        "p1_read_snps": None,
        "p2_read_snps": None,
        "ref_ld_snps": None,
        "weight_ld_snps": None,
        "merge_ref_snps": None,
        "merge_reg_snps": None,
        "merge_sumstats_snps": None,
        "valid_alleles_snps": None,

        "p1_h2_obs": None,
        "p1_h2_obs_se": None,
        "p1_lambda_gc": None,
        "p1_mean_chi2": None,
        "p1_intercept": None,
        "p1_intercept_se": None,
        "p1_ratio": None,
        "p1_ratio_se": None,

        "p2_h2_obs": None,
        "p2_h2_obs_se": None,
        "p2_lambda_gc": None,
        "p2_mean_chi2": None,
        "p2_intercept": None,
        "p2_intercept_se": None,
        "p2_ratio": None,
        "p2_ratio_se": None,

        "gencov": None,
        "gencov_se": None,
        "mean_z1z2": None,
        "gencov_intercept": None,
        "gencov_intercept_se": None,

        "rg": None,
        "rg_se": None,
        "rg_z": None,
        "rg_p": None,

        "summary_p1_path": None,
        "summary_p2_path": None,
        "summary_rg": None,
        "summary_rg_se": None,
        "summary_rg_z": None,
        "summary_rg_p": None,
        "summary_h2_obs": None,
        "summary_h2_obs_se": None,
        "summary_h2_int": None,
        "summary_h2_int_se": None,
        "summary_gcov_int": None,
        "summary_gcov_int_se": None,

        "has_analysis_end": False,
        "has_rg": False,
        "has_gencov": False,
        "has_p1_h2": False,
        "has_p2_h2": False,
        "has_summary_table": False,
        "finished_successfully": False,
    }

    row["analysis_start"] = find_first(RE_BEGINNING_ANALYSIS, text)
    row["analysis_end"] = find_first(RE_FINISHED_ANALYSIS, text)
    row["total_time_sec"] = safe_float(find_first(RE_TOTAL_TIME, text))

    row["has_analysis_end"] = row["analysis_end"] is not None

    row["out_prefix"] = find_first(RE_OUT, text)
    row["ref_ld_chr"] = find_first(RE_REF_LD, text)
    row["w_ld_chr"] = find_first(RE_W_LD, text)
    row["rg_arg_raw"] = find_first(RE_RG, text)

    if row["rg_arg_raw"]:
        rg_parts = [x.strip() for x in str(row["rg_arg_raw"]).split(",")]
        if len(rg_parts) >= 1:
            row["p1_path"] = rg_parts[0]
            row["p1_id"] = basename_without_sumstats(rg_parts[0])
        if len(rg_parts) >= 2:
            row["p2_path"] = rg_parts[1]
            row["p2_id"] = basename_without_sumstats(rg_parts[1])

    read_sumstats_files = find_all(RE_READ_SUMSTATS, text)
    row["summary_stats_files_n"] = len(read_sumstats_files)

    if read_sumstats_files:
        if row["p1_path"] is None and len(read_sumstats_files) >= 1:
            row["p1_path"] = read_sumstats_files[0]
            row["p1_id"] = basename_without_sumstats(read_sumstats_files[0])
        if row["p2_path"] is None and len(read_sumstats_files) >= 2:
            row["p2_path"] = read_sumstats_files[1]
            row["p2_id"] = basename_without_sumstats(read_sumstats_files[1])

    snp_reads = find_all(RE_READ_SUMSTATS_SNPS, text)
    if len(snp_reads) >= 1:
        row["p1_read_snps"] = safe_int(snp_reads[0])
    if len(snp_reads) >= 2:
        row["p2_read_snps"] = safe_int(snp_reads[1])

    row["ref_ld_snps"] = safe_int(find_first(RE_REF_SNPS, text))
    row["weight_ld_snps"] = safe_int(find_first(RE_WEIGHT_SNPS, text))
    row["merge_ref_snps"] = safe_int(find_first(RE_AFTER_MERGE_REF, text))
    row["merge_reg_snps"] = safe_int(find_first(RE_AFTER_MERGE_REG, text))
    row["merge_sumstats_snps"] = safe_int(find_first(RE_AFTER_MERGE_SUMSTATS, text))
    row["valid_alleles_snps"] = safe_int(find_first(RE_VALID_ALLELES, text))

    m_rg_comp = RE_COMPUTING_RG.search(text)
    if m_rg_comp:
        row["computing_rg_current"] = safe_int(m_rg_comp.group(1))
        row["computing_rg_total"] = safe_int(m_rg_comp.group(2))

    m1 = RE_PHENO1_H2.search(text)
    if m1:
        row["p1_h2_obs"] = safe_float(m1.group(1))
        row["p1_h2_obs_se"] = safe_float(m1.group(2))
        row["p1_lambda_gc"] = safe_float(m1.group(3))
        row["p1_mean_chi2"] = safe_float(m1.group(4))
        row["p1_intercept"] = safe_float(m1.group(5))
        row["p1_intercept_se"] = safe_float(m1.group(6))
        row["p1_ratio"] = safe_float(m1.group(7))
        row["p1_ratio_se"] = safe_float(m1.group(8))
        row["has_p1_h2"] = True

    m2 = RE_PHENO2_H2.search(text)
    if m2:
        row["p2_h2_obs"] = safe_float(m2.group(1))
        row["p2_h2_obs_se"] = safe_float(m2.group(2))
        row["p2_lambda_gc"] = safe_float(m2.group(3))
        row["p2_mean_chi2"] = safe_float(m2.group(4))
        row["p2_intercept"] = safe_float(m2.group(5))
        row["p2_intercept_se"] = safe_float(m2.group(6))
        row["p2_ratio"] = safe_float(m2.group(7))
        row["p2_ratio_se"] = safe_float(m2.group(8))
        row["has_p2_h2"] = True

    mg = RE_GENCOV.search(text)
    if mg:
        row["gencov"] = safe_float(mg.group(1))
        row["gencov_se"] = safe_float(mg.group(2))
        row["mean_z1z2"] = safe_float(mg.group(3))
        row["gencov_intercept"] = safe_float(mg.group(4))
        row["gencov_intercept_se"] = safe_float(mg.group(5))
        row["has_gencov"] = True

    mc = RE_GENCOR.search(text)
    if mc:
        row["rg"] = safe_float(mc.group(1))
        row["rg_se"] = safe_float(mc.group(2))
        row["rg_z"] = safe_float(mc.group(3))
        row["rg_p"] = safe_float(mc.group(4))
        row["has_rg"] = True

    ms = RE_SUMMARY_ROW.search(text)
    if ms:
        row["summary_p1_path"] = ms.group(1)
        row["summary_p2_path"] = ms.group(2)
        row["summary_rg"] = safe_float(ms.group(3))
        row["summary_rg_se"] = safe_float(ms.group(4))
        row["summary_rg_z"] = safe_float(ms.group(5))
        row["summary_rg_p"] = safe_float(ms.group(6))
        row["summary_h2_obs"] = safe_float(ms.group(7))
        row["summary_h2_obs_se"] = safe_float(ms.group(8))
        row["summary_h2_int"] = safe_float(ms.group(9))
        row["summary_h2_int_se"] = safe_float(ms.group(10))
        row["summary_gcov_int"] = safe_float(ms.group(11))
        row["summary_gcov_int_se"] = safe_float(ms.group(12))
        row["has_summary_table"] = True

    if row["has_analysis_end"] and row["has_rg"]:
        row["finished_successfully"] = True

    m_error = RE_ERROR.search(text)
    if m_error:
        row["status"] = "ERROR"
        row["error_message"] = m_error.group(1).strip().replace("\n", " | ")

    if not row["finished_successfully"] and row["status"] != "ERROR":
        row["status"] = "INCOMPLETE"

    return row


#####################
# file search logic #
#####################

def collect_log_files(
    input_dir: Path,
    recursive: bool = True,
    sort_files: bool = True
) -> List[Path]:
    if recursive:
        files = list(input_dir.rglob("*.log"))
    else:
        files = list(input_dir.glob("*.log"))

    if sort_files:
        files = sorted(files)

    return files


##############################
# save paths grouped by type #
##############################

def save_status_file_lists(rows: List[Dict[str, Optional[object]]], output_tsv: Path) -> None:
    base_path = output_tsv.with_suffix("")

    status_files = {
        "OK": base_path.with_name(base_path.name + "_OK.txt"),
        "ERROR": base_path.with_name(base_path.name + "_ERROR.txt"),
        "INCOMPLETE": base_path.with_name(base_path.name + "_INCOMPLETE.txt"),
        "PARSER_ERROR": base_path.with_name(base_path.name + "_PARSER_ERROR.txt"),
    }

    handles = {}
    try:
        for status, path in status_files.items():
            handles[status] = open(path, "w", encoding="utf-8")

        for row in rows:
            status = row.get("status", "UNKNOWN")
            log_file = row.get("log_file")

            if status in handles and log_file is not None:
                handles[status].write(f"{log_file}\n")
    finally:
        for handle in handles.values():
            handle.close()

    print("\nSaved file lists:")
    for status, path in status_files.items():
        print(f"  {status}: {path}")


#################
# main function #
#################

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse LDSC rg log files into one TSV table."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing LDSC .log files"
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="Output TSV file"
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search recursively for .log files"
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Maximum number of first log files to process (for testing). Default: all files."
    )
    parser.add_argument(
        "--no-sort-files",
        action="store_true",
        help="Do not sort file list before processing."
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_tsv = Path(args.output_tsv)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    log_files = collect_log_files(
        input_dir=input_dir,
        recursive=args.recursive,
        sort_files=not args.no_sort_files
    )

    if not log_files:
        raise FileNotFoundError(f"No .log files found in: {input_dir}")

    if args.max_files is not None:
        if args.max_files <= 0:
            raise ValueError("--max-files must be > 0")
        log_files = log_files[:args.max_files]

    print(f"Files to parse: {len(log_files)}")

    rows: List[Dict[str, Optional[object]]] = []

    for log_file in tqdm(log_files, desc="Parsing LDSC logs", unit="log"):
        try:
            row = parse_ldsc_log(log_file)
        except Exception as e:
            row = {
                "log_file": str(log_file),
                "log_filename": log_file.name,
                "status": "PARSER_ERROR",
                "error_message": str(e),

                "analysis_start": None,
                "analysis_end": None,
                "total_time_sec": None,

                "out_prefix": None,
                "ref_ld_chr": None,
                "w_ld_chr": None,
                "rg_arg_raw": None,

                "p1_path": None,
                "p2_path": None,
                "p1_id": None,
                "p2_id": None,

                "summary_stats_files_n": None,
                "computing_rg_current": None,
                "computing_rg_total": None,

                "p1_read_snps": None,
                "p2_read_snps": None,
                "ref_ld_snps": None,
                "weight_ld_snps": None,
                "merge_ref_snps": None,
                "merge_reg_snps": None,
                "merge_sumstats_snps": None,
                "valid_alleles_snps": None,

                "p1_h2_obs": None,
                "p1_h2_obs_se": None,
                "p1_lambda_gc": None,
                "p1_mean_chi2": None,
                "p1_intercept": None,
                "p1_intercept_se": None,
                "p1_ratio": None,
                "p1_ratio_se": None,

                "p2_h2_obs": None,
                "p2_h2_obs_se": None,
                "p2_lambda_gc": None,
                "p2_mean_chi2": None,
                "p2_intercept": None,
                "p2_intercept_se": None,
                "p2_ratio": None,
                "p2_ratio_se": None,

                "gencov": None,
                "gencov_se": None,
                "mean_z1z2": None,
                "gencov_intercept": None,
                "gencov_intercept_se": None,

                "rg": None,
                "rg_se": None,
                "rg_z": None,
                "rg_p": None,

                "summary_p1_path": None,
                "summary_p2_path": None,
                "summary_rg": None,
                "summary_rg_se": None,
                "summary_rg_z": None,
                "summary_rg_p": None,
                "summary_h2_obs": None,
                "summary_h2_obs_se": None,
                "summary_h2_int": None,
                "summary_h2_int_se": None,
                "summary_gcov_int": None,
                "summary_gcov_int_se": None,

                "has_analysis_end": False,
                "has_rg": False,
                "has_gencov": False,
                "has_p1_h2": False,
                "has_p2_h2": False,
                "has_summary_table": False,
                "finished_successfully": False,
            }

        rows.append(row)

    all_columns: List[str] = []
    seen = set()

    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                all_columns.append(key)

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    with output_tsv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=all_columns,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved parsed results to: {output_tsv}")

    status_counts = {}
    for row in rows:
        status = row.get("status", "UNKNOWN")
        status_counts[status] = status_counts.get(status, 0) + 1

    print("Status summary:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")

    save_status_file_lists(rows=rows, output_tsv=output_tsv)


if __name__ == "__main__":
    main()
