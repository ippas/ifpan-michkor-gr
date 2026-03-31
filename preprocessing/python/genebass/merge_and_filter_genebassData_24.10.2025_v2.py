import os
import glob
import time
import psutil
import gzip
import polars as pl
from tqdm import tqdm


def filter_polars_table(df: pl.DataFrame, verbose: bool = True, **filters) -> pl.DataFrame:
    """
    Filter a Polars DataFrame using flexible keyword-based conditions.
    Example:
        filter_polars_table(df, pvalue__lt=0.0001, annotation__eq="pLoF")
    """

    if not filters:
        if verbose:
            print("⚠️ No filters provided — returning original DataFrame.")
        return df

    # Build expressions
    exprs = []
    for key, val in filters.items():
        if "__" not in key:
            raise ValueError(f"Invalid filter '{key}'. Use syntax like 'pvalue__lt=0.05'")

        col, op = key.split("__", 1)
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in DataFrame.")

        c = pl.col(col)
        match op:
            case "lt": exprs.append(c < val)
            case "le": exprs.append(c <= val)
            case "gt": exprs.append(c > val)
            case "ge": exprs.append(c >= val)
            case "eq": exprs.append(c == val)
            case "ne": exprs.append(c != val)
            case _: raise ValueError(f"Unsupported operator: '{op}'")

    condition = exprs[0]
    for e in exprs[1:]:
        condition &= e

    before_rows = df.height
    filtered = df.filter(condition)
    after_rows = filtered.height

    if verbose:
        reduction = 100 * (1 - after_rows / before_rows) if before_rows > 0 else 0
        print(f"✅ Filtered: {after_rows:,} rows (from {before_rows:,}) → {reduction:.1f}% reduction")

    return filtered


def merge_filtered_bgz_files(
    input_directory: str,
    output_path: str,
    file_to_exclude: str | None = None,
    memory_limit_gb: float = 80.0,
    max_files: int | None = None,
    return_polars_df: bool = False,
    verbose: bool = True,
    **filters
) -> pl.DataFrame | None:
    """
    Merge and filter multiple .tsv.bgz files into a single gzip-compressed output file.
    Optionally returns the merged Polars DataFrame.

    Parameters
    ----------
    input_directory : str
        Directory with .tsv.bgz files.
    output_path : str
        Path to save the merged filtered file (.bgz).
    file_to_exclude : str | None
        File to skip (e.g., already merged file).
    memory_limit_gb : float
        Stop the process if memory usage exceeds this limit.
    max_files : int | None
        Limit to first N files for testing/debugging.
    return_polars_df : bool, default False
        If True, return merged Polars DataFrame.
    verbose : bool
        Print progress and stats.
    **filters :
        Filters like pvalue__lt=0.0001 or annotation__eq="pLoF".
    """

    # --- 1️⃣ Find files ---
    input_files = sorted(glob.glob(os.path.join(input_directory, "*.tsv.bgz")))
    if file_to_exclude:
        input_files = [f for f in input_files if not f.endswith(file_to_exclude)]

    if not input_files:
        raise FileNotFoundError(f"No .tsv.bgz files found in {input_directory}")

    if max_files is not None:
        input_files = input_files[:max_files]
        print(f"🧪 Test mode: limiting to first {len(input_files)} files.")

    if verbose:
        print(f"📂 Found {len(input_files)} files to merge.")
        if file_to_exclude:
            print(f"🚫 Excluding file: {file_to_exclude}")

    start_time = time.time()
    partial_results = []
    total_rows_before = 0
    total_rows_after = 0

    # --- 2️⃣ Process each file ---
    for file in tqdm(input_files, desc="🔄 Filtering files", unit="file"):
        mem_gb = psutil.virtual_memory().used / (1024 ** 3)
        if mem_gb > memory_limit_gb:
            raise MemoryError(f"❌ Memory limit exceeded ({mem_gb:.2f} GB > {memory_limit_gb} GB).")

        try:
            df = pl.read_csv(
                file,
                separator="\t",
                has_header=True,
                infer_schema_length=100,
                null_values=["NA", "null", "NaN"],
                low_memory=True
            )
        except Exception as e:
            print(f"⚠️ Skipping file {file} due to read error: {e}")
            continue

        total_rows_before += df.height

        if filters:
            df = filter_polars_table(df, verbose=False, **filters)

        total_rows_after += df.height

        if df.height == 0:
            continue

        partial_results.append(df)

        if verbose:
            print(f"✅ {os.path.basename(file)} → {df.height:,} rows kept")

        mem_gb = psutil.virtual_memory().used / (1024 ** 3)
        if mem_gb > memory_limit_gb:
            raise MemoryError(f"❌ Memory limit exceeded ({mem_gb:.2f} GB > {memory_limit_gb} GB).")

    if not partial_results:
        print("⚠️ No data left after filtering — nothing to save.")
        return None

    merged_df = pl.concat(partial_results, how="vertical_relaxed")
    del partial_results

    # --- 6️⃣ Save final merged file (.bgz) ---
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    csv_content = merged_df.write_csv(None, separator="\t", include_header=True)
    csv_text = csv_content.decode("utf-8") if isinstance(csv_content, bytes) else csv_content

    with gzip.open(output_path, "wt", encoding="utf-8") as fout:
        fout.write(csv_text)

    elapsed = time.time() - start_time
    print("\n✅ Merge & filtering completed.")
    print(f"📄 Output saved to: {output_path}")
    print(f"📊 Rows: {total_rows_before:,} → {total_rows_after:,} after filtering")
    print(f"⏱️ Time: {elapsed/60:.2f} minutes")

    if return_polars_df:
        return merged_df
    else:
        del merged_df
        return None



df = merge_filtered_bgz_files(
    input_directory="/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/all_categories_SKATO",
    output_path="/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/merged_and_filtered_genebassData/genebass_allCategories_SKATO_mergedFilteredP01.tsv.bgz",
    file_to_exclude="genebass_merged_allCategories_SKATO.tsv.bgz",
    memory_limit_gb=80,
    pvalue__lt=0.01,
    max_files=None,
    return_polars_df=True,
    verbose=False
)
