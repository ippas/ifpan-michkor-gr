import hail as hl
import pandas as pd
import numpy as np
from pyspark.sql import SparkSession
import os

from uuid import uuid4
import platform
import socket
from pyspark.sql import SparkSession

hl_init_kwargs = {
    'default_reference': 'GRCh38'
}

def hl_init(**kwargs):
    hl_init_kwargs.update(kwargs)
    hl.init(**hl_init_kwargs)

hl_init_kwargs = {
    'default_reference': 'GRCh38'
}


hl_init()
mt = hl.read_matrix_table('/net/pr2/projects/plgrid/plggneuromol/resources/genebass-500k/results.mt')

filtered_mt = mt.filter_rows(
    (mt.gene_symbol == "TSPAN6")
)

filtered_mt.count()

filtered_mt.show()

filtered_mt.describe()

filtered_mt.describe().show()


mt.describe()

mt.phenocode.describe()


mt.phenocode.show()

mt.phenocode.show(100)

mt.rows().select().show()

mt.phenocode.show(100)



mt.total_variants_pheno.describe()

mt.describe()

mt.cols().show()
mt.rows().show()
mt.col.describe()


mt.describe()
filtered_mt = mt.filter_cols(mt.phenocode == "C15")
filtered_mt.count()


# Select distinct gene_symbols from the filtered MatrixTable
gene_symbols = filtered_mt.rows().select(filtered_mt.gene_symbol).distinct().collect()

len(filtered_mt.gene_symbol.collect())

len(filtered_mt.Pvalue.collect())

# Filter rows where any Pvalue is less than 0.01
significant_results = filtered_mt.filter_rows(
    hl.agg.any(filtered_mt.Pvalue < 0.00001)
)

significant_results.count_rows()

significant_results.gene_symbol.show()
