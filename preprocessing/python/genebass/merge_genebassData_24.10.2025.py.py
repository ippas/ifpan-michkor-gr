import os
import gzip
import glob
import time
import psutil
from tqdm import tqdm

def merge_bgz_files(
    input_directory: str,
    output_path: str,
    file_to_exclude: str | None = None,
    memory_limit_gb: float = 100.0
):
    """
    Merge multiple .tsv.bgz files from a given directory into a single output file.
    The function supports progress tracking, runtime measurement, 
    and memory usage monitoring with an optional limit.

    Parameters
    ----------
    input_directory : str
        Path to the directory containing .tsv.bgz files.
    output_path : str
        Full path to the output file (.tsv.bgz).
    file_to_exclude : str | None, default None
        File name to be excluded from merging (e.g., 'merged_all.tsv.bgz').
    memory_limit_gb : float, default 100.0
        Maximum allowed memory usage (in GB). Process will stop if exceeded.
    """

    # 🔍 Find all .tsv.bgz files in directory
    input_files = sorted(glob.glob(os.path.join(input_directory, "*.tsv.bgz")))
    if file_to_exclude:
        input_files = [f for f in input_files if not f.endswith(file_to_exclude)]

    if not input_files:
        raise FileNotFoundError(f"No .tsv.bgz files found in {input_directory}")

    print(f"📂 Found {len(input_files)} files to merge.")
    print(f"🚫 Excluded file: {file_to_exclude if file_to_exclude else 'None'}")

    start_time = time.time()
    header_written = False

    # ✍️ Open the output file for writing (gzip-compressed)
    with gzip.open(output_path, "wt") as fout:
        for i, file in enumerate(tqdm(input_files, desc="🔄 Merging files", unit="file")):
            
            # 🔒 Check memory usage
            mem_gb = psutil.virtual_memory().used / (1024 ** 3)
            if mem_gb > memory_limit_gb:
                raise MemoryError(
                    f"❌ Memory limit exceeded ({mem_gb:.2f} GB > {memory_limit_gb} GB). Process aborted."
                )

            with gzip.open(file, "rt") as fin:
                for j, line in enumerate(fin):
                    # Write the header only once
                    if not header_written:
                        fout.write(line)
                        header_written = True
                    elif j > 0:  # Skip headers from subsequent files
                        fout.write(line)

    elapsed = time.time() - start_time
    print("\n✅ Merge completed successfully.")
    print(f"📄 Output saved to: {output_path}")
    print(f"⏱️ Total time: {elapsed/60:.2f} minutes ({elapsed:.1f} seconds)")
    print(f"📊 Files merged: {len(input_files)}")


merge_bgz_files(
    input_directory="/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/all_categories_SKATO",
    output_path="/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/all_categories_SKATO/genebass_merged_allCategories_SKATO_24.10.2025.tsv.bgz",
    memory_limit_gb=100
)

merge_bgz_files(
    input_directory="/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/all_categories_SKAT",
    output_path="/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/all_categories_SKAT/genebass_merged_allCategories_SKAT_24.10.2025.tsv.bgz",
    memory_limit_gb=100
)

merge_bgz_files(
    input_directory="/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/all_categories_burden",
    output_path="/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/all_categories_burden/genebass_merged_allCategories_burden_24.10.2025.tsv.bgz",
    memory_limit_gb=100
)