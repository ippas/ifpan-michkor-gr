#!/bin/bash
###############################################################
### SLURM preamble:  ###
###############################################################
#SBATCH --account plgdepresja2-cpu
#SBATCH --partition plgrid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time 72:00:00
#SBATCH -C localfs
#############################################################

# Extraction of following demographic data from UKb phenotypes:
# column code in html - column code in txt  - Description
# 22  31  Sex
# 23  34  Year of birth
# 6488-6489 12144 Height
# 9751-9754 21002 Weight
# 9747-9750 21001 BMI
# 72-75 48  Waist circumference
# 76-79 49  Hip circumference
# 1550-1557 4079  Diastolic blood pressure
# 1558-1565 4080  Systolic blood pressure
# 657-660 1259  Smoking/smokers in household ***************************************************** 
# 661-664 1269  Exposure to tobacco smoke at home  ***********************************************
# 665-668 1279  Exposure to tobacco smoke outside home  ******************************************
# 1076-1079 2644 Light smokers  ******************************************************************
# 1343-1346 3159  Smoked cigarette or pipe within last hour  *************************************
# 4388-4391 5959  Previously smoked cigarettes on most/all days  *********************************
# 769-772 1558  Alcohol intake frequency
# 793-796 1618  Alcohol usually taken with meals
# 1450-1453 3731  Former alcohol drinker  ********************************************************
# 1048-1051 2443  Diabetes diagnosed by doctor  **************************************************
# 1196-1199 2976  Age diabetes diagnosed  ********************************************************
# 1538-1541 4041  Gestational diabetes only ******************************************************
# 6120-6135 6153  Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones  *************
# 6136-6159 6154  Medication for pain relief, constipation, heartburn  ***************************
# 6160-6187 6155  Vitamin and mineral supplements  ***********************************************
# 6188-6203 6156  Manic/hyper symptoms  **********************************************************
# 6320-6331 6177  Medication for cholesterol, blood pressure or diabetes  ************************
# 452-455 137 Number of treatments/medications taken
# 1064-1067 2492  Taking other prescription medications  *****************************************
# 9623  20466 Ever prescribed a medication for unusual or psychotic experiences ******************
# 9729-9731 20551 Substance of prescription or over-the-counter medication addiction  ************
# 1052-1055 2453  Cancer diagnosed by doctor *****************************************************
# 14656-14677 40008 Age at cancer diagnosis ******************************************************
# 13802-13804 30080 Platelet count
# 13867-13869 30130 Monocyte count
# 13711-13713 30010 Red blood cell (erythrocyte) count
# 13698-13700 30000 White blood cell (leukocyte) count
# 13750-13752 30040 Mean corpuscular volume
# 13763-13765 30050	Mean corpuscular haemoglobin
# 13828-13830 30100 Mean platelet (thrombocyte) volume
# 93-96 54  UK Biobank assessment centre
# 8908  20118 Home area population density - urban or rural

cat output_ukb.txt | sed -e 's/^\t*//' | cut -f 1,23,24,6489-6490,9752-9755,9748-9751,73-76,77-80,1551-1558,1559-1566,658-661,662-665,666-669,1077-1080,1344-1347,4389-4392,770-773,794-797,1451-1454,1049-1052,1197-1200,1539-1542,6121-6136,6137-6160,6161-6188,6189-6204,6321-6332,453-456,1065-1068,9624,9730-9732,1053-1056,14657-14678,13803-13805,13868-13870,13712-13714,13699-13701,13751-13753,13764-13766,13829-13831,94-97,8909 > final.tsv

exit;
