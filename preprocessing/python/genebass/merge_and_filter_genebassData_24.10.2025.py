import gzip
import polars as pl

def load_bgz_polars(file_path: str) -> pl.DataFrame:
    import polars as pl
    print(f"📥 Loading: {file_path}")
    df = pl.read_csv(
        file_path,
        separator="\t",
        has_header=True,
        infer_schema_length=100,
        null_values=["NA", "null", "NaN"],
        low_memory=True
    )
    print(f"✅ {df.shape[0]:,} rows × {df.shape[1]} columns")
    return df

import polars as pl

def filter_polars_table(
    df: pl.DataFrame,
    verbose: bool = True,
    **filters
) -> pl.DataFrame:
    """
    Filter a Polars DataFrame using flexible keyword-based conditions.
    Example usage:
        filter_polars_table(df, pvalue__lt=0.0001, annotation__eq="pLoF")

    Supported operators:
        __lt  → less than (<)
        __le  → less or equal (≤)
        __gt  → greater than (>)
        __ge  → greater or equal (≥)
        __eq  → equal (==)
        __ne  → not equal (!=)

    Parameters
    ----------
    df : pl.DataFrame
        The input Polars DataFrame.
    verbose : bool, default True
        If True, print filtering summary.
    **filters :
        Filtering conditions as keyword arguments.

    Returns
    -------
    pl.DataFrame
        The filtered Polars DataFrame.
    """

    if not filters:
        if verbose:
            print("⚠️ No filters provided — returning original DataFrame.")
        return df

    # Build expressions for each filter
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

    # Combine all expressions with logical AND
    condition = exprs[0]
    for e in exprs[1:]:
        condition &= e

    before_rows = df.height
    before_size = df.estimated_size()

    # Apply filter
    filtered = df.filter(condition)

    after_rows = filtered.height
    after_size = filtered.estimated_size()

    if verbose:
        reduction = 100 * (1 - after_rows / before_rows) if before_rows > 0 else 0
        size_diff_mb = (before_size - after_size) / (1024**2)
        print(f"✅ Filtered: {after_rows:,} rows (from {before_rows:,}) → {reduction:.1f}% reduction")
        print(f"📉 Estimated memory saved: {size_diff_mb:.2f} MB")

    return filtered



df = load_bgz_polars("/home/mateusz/projects/ifpan-michkor-gr/data/genebass_v2/all_categories_SKATO/genebass_SKATO_UK_Biobank_Assessment_Centre_-_Verbal_interview_-_Operations.tsv.bgz")

filtered = filter_polars_table(df, pvalue__lt=0.0001)

