#!/usr/bin/env bash

set -euo pipefail

BASE_DIR="/home/mateusz/projects/ifpan-michkor-gr"
DATA_DIR="${BASE_DIR}/data/dbsnp_rsidAnnotation/raw"

mkdir -p "${DATA_DIR}"
cd "${DATA_DIR}"

# full dbSNP for GRCh37
wget -c https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget -c https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi

# optional: smaller file with common variants only
# wget -c https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
# wget -c https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz.tbi

echo "dbSNP files are ready in: ${DATA_DIR}"
