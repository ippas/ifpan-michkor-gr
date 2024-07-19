#%%
#SBATCH --partition plgrid-now
#SBATCH --ntasks-per-node=2
#SBATCH --mem=110GB
#SBATCH --time 12:00:00
#SBATCH --job-name jupyter-subjob

from pathlib import Path
from math import ceil

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

WD = Path(
    '/net/pr2/projects/plgrid/plggneuromol/Alireza', 
    'UkbPhenotypes/All_required_data/Binary_Standardization',
    'PCA_Data/New_PCAs/PWRSM/'
)


#%% read covariate data and genetic data columns
cov_data_path = WD / 'merged_covariate_table.tsv'
gen_data_path = WD / 'UKBgeneticPCandACGWASEuropeans_filtered_without_last.tsv'

cov_data = pd.read_csv(cov_data_path, sep='\t')
cov_cols = [
    '2200901', '2200902', '2200903', '2200904', '2200905', '2200906',
    '2200907', '2200908', '2200909', '22009010', '22009011', '22009012',
    '22009013', '22009014', '22009015', '22009016', '22009017', '22009018',
    '22009019', '22009020'
]
cov_data = cov_data[cov_cols]

gen_data_cols = pd.read_csv(gen_data_path, sep='\s+', nrows=0)

#%% adjust genetic data
X = cov_data.to_numpy()
model = LinearRegression()

max_cols = 1_702
n =  gen_data_cols.shape[1]
part_paths = []
for i in range(ceil(n / max_cols)):
    print(f'Part {i}')
    gen_slice = range(i * max_cols, min(n, (i + 1) * max_cols))

    y = np.genfromtxt(
        gen_data_path, dtype=float, usecols=gen_slice, skip_header=1
    )
    y_pred = model.fit(X, y).predict(X)
    y_adj = y - y_pred

    part_path = f'adjusted-matrix-data-{gen_slice.start}-{gen_slice.stop}.tsv'
    np.savetxt(
        part_path,
        y_adj,
        delimiter='\t',
        fmt='%.9f',
        header=gen_data_cols.iloc[:, gen_slice].to_csv(sep='\t', index=False).strip(),
        comments='',
    )
    part_paths.append(part_path)


#%% merge files
print('Merging files')
output_file = 'adjusted_genetic_data.tsv'
files = [open(path, 'r') for path in part_paths]

with open(output_file, 'w') as out:
    for lines in zip(*files):
        merged_line = '\t'.join(line.strip() for line in lines) + '\n'
        out.write(merged_line)
    for file in files:
        file.close()
