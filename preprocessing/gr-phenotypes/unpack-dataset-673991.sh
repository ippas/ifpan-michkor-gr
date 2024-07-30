#!/usr/bin/env bash

#SBATCH -A plgwgsdepresja3-cpu
#SBATCH --partition plgrid
#SBATCH --mem=20GB
#SBATCH --time 1-00:00:00
#SBATCH --job-name ukb-r
#SBATCH --output=/net/pr2/projects/plgrid/plggneuromol/matzieb/projects/ifpan-michkor-gr/tools/slurm-log/%j.out
#SBATCH --error=/net/pr2/projects/plgrid/plggneuromol/matzieb/projects/ifpan-michkor-gr/tools/slurm-log/%j.err


ID="673991"

echo "MD5..."
tools/ukb-utils/ukbmd5 raw/ukb-phenotypes/ukb$ID.enc

echo "Decrypt..."
tools/ukb-utils/ukbunpack raw/ukb-phenotypes/ukb$ID.enc raw/ukb-phenotypes/k62979r$ID.key

echo "Convert..."
tools/ukb-utils/ukbconv raw/ukb-phenotypes/ukb$ID.enc_ukb docs -etools/ukb-utils/encoding.ukb
tools/ukb-utils/ukbconv raw/ukb-phenotypes/ukb$ID.enc_ukb r -etools/ukb-utils/encoding.ukb