## UKBiobank Phenotypic Data Processing

* Almost all codes and scripts in this document ran on HPC system with slurm job scheduler through following cinfiguration:

```
#SBATCH --nodes=1
#SBATCH --mem=180G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time [various]
```

### Stage1 (Raw Data Processing)

After downloading **packed UKB archive file**, the main data extracted through the utilization of ***[ukbunpack](https://biobank.ndph.ox.ac.uk/showcase/util/ukbunpack).***

Next, ***[ukbconv](https://biobank.ndph.ox.ac.uk/showcase/util/ukbconv)*** is used for converting the raw data to desired format like csv, txt, html, etc.

> **html** file must be extracted in anyway as it has all the avilable data in a well-organized format for reviewing and data parsing.

After parsing UKB data and obtaining **output_ukb.txt** (which contains all the information from original phenotypic data) [approximately 18.5GB], the **html** file has been checked and looked for the examples of ICD10 codes. Then, cut command applied to find the columns in the original data which contains **an example of ICD10** codes.

```
cat output_ukb.txt | cut -f 14800-15000 | grep X998 | head
```

After finding the interested columns’ number, the command line for extraction of **JUST ICD10** codes (not ICD10 codes mixed with other words) for each patient would be the following:

```
cat output_ukb.txt | cut -f 2,15811-16069 > RealDiseaseCodesForPatients.txt
```

**RealDiseaseCodesForPatients.txt** contains all the disease codes for each patient, indicating what participant has what disorder. on the other side, the file **ICD_10s.tsv** (all ICD10 names extracted from **html** file) was used in python code *1_UniqueICDs.py* to have **unique_disease_codes.tsv** to remove duplicates from the list.

> The next step was a crucial data processing practice for preparing a **matrix** of samples and diseasecodes with 1 for the presence of the code and 0 for the absent of the code.

We removed column names (which are just numbers) from the result of above code and just have the sample IDs with disease codes for each of them. This file would be suitable for making the desired **matrix**. The related code for making the **matrix** is brought in ***2_CreatingMatrixFromRealDiseaseCodes.py.***

At the same time, the frequency of each disease code in UKB participants calculated by ***3_ICD10RealFrequencies.py*** to know about the number of samples who had the specific disorder. The result used for quick check for any discrepancy between the frequencies in **unique_disease_code_frequencies.txt** from ***3_ICD10RealFrequencies.py*** and the created matrix, named as **output_file4.csv**.

The first column in **ouput_file4.csv** contianed sample IDs and it was removed as it could cause some dificulties in data normalization step.

```
cut -d',' -f2- output_file4.csv > output_file4_withoutEID.csv

```

Since, the data was too large for doing the normalization, we splited that to "500samples" chunks and prepared them for the standardization step by ***5_SplitCreatedMatrixForLaterNormalization.sh***.

Data normalization performed on all chunks separately (in parallel), using ***6_DataNormalization.py*** and then all the standardized data merged into a one file and extra empty columns removed by ***7_MergingNormalizedICD10Data.sh***

> The output file was ready to add to the common GWAS phenotypic data for each samples.

Selected common GWAS (demographic) phenotypes and related codes in UKB data for the participants in current study included:

```
31 Sex
34 Year of birth
48 Waist circumference
49 Hip circumference
54 UK Biobank assessment centre
137 Number of treatments/medications taken
1259 Smoking/smokers in household
1558 Alcohol intake frequency
1618 Alcohol usually taken with meals
2443 Diabetes diagnosed by doctor
2453 Cancer diagnosed by doctor
2492 Taking other prescription medications
4041 Gestational diabetes only
4079 Diastolic blood pressure
4080 Systolic blood pressure
6153 Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones
6154 Medication for pain relief, constipation, heartburn
6155 Vitamin and mineral supplements
6156 Manic/hyper symptoms
6177 Medication for cholesterol, blood pressure, or diabetes
12144 Height
20466 Ever prescribed a medication for unusual or psychotic experiences
20551 Substance of prescription or over-the-counter medication addiction
21001 BMI
21002 Weight
30000 White blood cell (leukocyte) count
30010 Red blood cell (erythrocyte) count
30040 Mean corpuscular volume
30050 Mean corpuscular haemoglobin
30080 Platelet count
30100 Mean platelet (thrombocyte) volume
30130 Monocyte count
40008 Age at cancer diagnosis
```

We first extracted the data by ***SlurmBashJobsForDemographic.sh***. However, some of thses phenotypes were categorized in more than one column as the number of times they did the assessment were two, three, four, etc.

In order to have just one unique column for each GWAS phenotype data, we used ***one-hot encoded*** approach:

- First, the categories with more than one column mergerd into array, using ***8_0_MakeMergeIntoArray.py*** for each category like waist circumference, vitamin and mineral supplement, etc.
- Next, ***one-hot encoded*** was done for each category separately by ***8_1_MakeOneHotEncoded.py***.
- Finally, the average data for the categories with more than one column calculated by the ***8_2_******GettingAverageForOneHotEncoded.py.***

> All the data have been merged and normalized into **AllGWASphenotypesNormalized.csv** and prepared to be added to our main **matrix**: **output_file4_withoutEID.csv**.

Before merging the main data, the selected samples from **AllGWASphenotypesNormalized.csv** must have been found and extracted from **output_file4.csv** as well. These samples were not considered for further studies as most of their GWAS phenotypes (from our list) had no available data! This way, all the related samples in our **matrix** would be matched with the data from GWAS phenotypes data. The utilized scripts could be found in **9_ExtractionOfRowNumbersFromICD10NormalizedMatchWithGWASsamples.py** and **10_ExtractInterestedRowsFromNormalizedICD.py**.

At the end of stage1, **AllGWASphenotypesNormalized.csv** and the output file from ***10_ExtractInterestedRowsFromNormalizedICD.py*** named as **FINALofFINALsICD10Normalized77_2.csv**.

> The output file (**ReadyForPCA.csv**) was ready to run PCA (principle component analysis).

### Stage 2 (Running PCA)

For the PCA, only the disease codes with frequency > 50 has been chosen from **unique_disease_code_frequencies.tsv** (converted from .txt to .tsv) by ***1_ExtractFreqMoreThan50.py***. The output file named as **filtered_statistics50.tsv**.

The result of statistics then utilized for filtering the main data (**ReadyForPCA.csv**) to have only the disease codes with frequency more than fifty.

Because of the data with relatively large size, a special tool named as ***[PCAone](https://github.com/Zilong-Li/PCAone)*** was employed for PCA calculations. Based on developer, ***PCAone*** is a fast and memory efficient PCA tool implemented in C++ aiming at providing comprehensive features and algorithms for different scenarios.

> We noticed the tool will transpose the input matrix based on sample size and features. Hence, we first transposed the main data and used it as input file for **PCAone**. In this case, the data will return to it's original structure after getting hit by **PCAone**. However, the transposing tool ***[Transpose](https://anaconda.org/molinerislab/transpose)*** we used, works with tsv file more accurately. So, the main data first converted from .csv to .tsv and after transposing, it returned to csv, ready to be as an input file to **PCAone**. This series of action performed by ***3_MatrixTransposing3.sh***.

PCA must be run on the numeric data with no headers as the combined words and numbers like what was seen in disease codes may interrupt the calculation. After preparing the input file for PCAone, first column (EID) and first row (disease codes) have been removed, using ***4_TsvToCsv.sh***.

**PCAone** accept the input file, compressed as [***ZSTD***](https://anaconda.org/conda-forge/zstd). The related code is provided in ***5_ZSTD.sh***.

PCAone performed (***6_PCAone.sh***) and the screeplot illustrated subsequently (***7_ScreePlot.py***).

Top phenotypes for first 10 components have been listed by pca.loadings and ***8_LoadingsAndTotalPhenotypes.py***. Then, first top 20 phenotypes from the extended top phenotypes list marked and extracted through ***9_MakingZlist.py***.

Next step included a series of data manipulation to provide a file cotaining all the GWAS and ICD10 codes with frequency > 50 as a single line header (***9_PhenotypeCodesHeader.sh***). We needed this file for creating a dictionary mapping codes from our top 20 phenotypes codes to code descriptions. The related script stored in ***10_FINALofFINALimportantPHENOTYPES2.py***.

Phenotypes with highest importance selected from each component and the related data for participants obtained form the transposed main matrix by ***11_NewLines.py*** and used for making PCA scatter plot as the final step for stage 2 (***12_PCAplotColored22.py***).

> All the above steps in stage 2 also performed for disease codes with frequency > 500 in order to have a comparison between the PCA results and identification of any potential new contributed phenotype(s).

### Stage 3 (Data Clustering)

In progress...
