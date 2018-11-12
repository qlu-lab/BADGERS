# BADGERS

BADGERS(Biobank-wide Association Discovery using GEnetic Risk Scores), a powerful method to perform polygenic score-based biobank-wide association scan.

## Prerequisites

The software is developed and tested in Linux and Mac OS environments.
- Python 2.7
- numpy
- scipy
- pandas

## Quick Start 

This section demonstrates the fundamental usage of BADGERS, which is to conduct the association between selected UK-biobank traits (total of 1738) and a complex disease sum stats. If user wants to run BADGERS with other complex traits, he needs to create corresponding [weight database](https://github.com/qlu-lab/BADGERS/wiki/Create-db-files) and [covariance](https://github.com/qlu-lab/BADGERS/wiki/Create-covariance-file) to replace UK-biobank traits as input in the sample below. 
### Step1: Downloads BADGERS

```
$ git clone https://github.com/qlu-lab/BADGERS
$ cd ./BADGERS
```
In order to test whether BADGERS function correctly in user's machine user can run [sample test](https://github.com/qlu-lab/BADGERS/wiki/Sample-test) to test BADGERS' function

### Step2: Downloads UK-biobank traits weight database and covariance
```
$ wget ftp://ftp.biostat.wisc.edu/pub/lu_group/BADGERS/UK_biobank_Round2.tar.gz
$ tar xf UK_biobank_Round2.tar.gz
```
This folder will include the following files/folders:
```
cov/      ## covariance file for traits
weight_db/      ## weight database for traits
UKbiobank_1738_inputlist.csv      ## name for traits
GWAS.txt  ## sample GWAS summary stats for testing purpose
```
UKbiobank_1738_inputlist.csv file contains ID and name for UK-Biobank traits

For each UK-Biobank trait there will be one corresponding database file in the weight_db folder and two corresponding covariance files in cov folder

For example for a trait with name 50_raw, there exist one database file 50_raw.db in weight_db folder that holds GWAS summary stats for this trait and two covariance files with name 50_raw_cov.txt.gz which held variance of SNPs in 50_raw.db and 50_raw_cov_tn.txt.gz which held phenotype of these SNPs in 1000 Genomes Project in cov folder


### Step3: Perform analysis

Here is how the user gets the association between UK_biobank traits and complex disease

```
python BADGERS.py \
--model_db_path UK_biobank_Round2/weight_db \
--covariance UK_biobank_Round2/cov \
--gwas_path UK_biobank_Round2/GWAS.txt \
--snp_column MarkerName \
--effect_allele_column Allele1 \
--non_effect_allele_column Allele2 \
--pvalue_column P.value \
--beta_column Effect \
--output_file output.csv \
--db_name 1738_traits.csv 
```
where
- *--model_db_path:*

    The path to folder contains database files for all input traits.

- *--covariance:*

    The path to folder contains covariance files for all input traits.

- *--gwas_path:*

    The location of input GWAS data(user input).

- *--snp_column:*

    Name of column holding SNP information.

- *--effect_allele_column:*

    Name of the column holding effect allele value.

- *--non_effect_allele_column:*

    Name of the column holding other/non effect allele value.

- *--pvalue_column:*

    Name of the column holding p-value.

- *--beta_column:*

    Name of the column holding beta value.

- *output_file:*

    Location where results will be saved.

- *db_name:*

    ID and name of selected traits involved in the analysis. To perform the analysis with all 1738 traits user can use file 1738_traits.csv as input. If user wants to perform analysis with specific traits user can state traits they want in the input file. For example, if user wants to test the association between 50_raw.db, 20016_raw.db and input disease, for db_name input user should have a selected_traits.csv file like:
<pre>
  50_raw,    Standing height
  20016_raw, Fluid intelligence score
</pre>

**Reminder**: In order to perform analysis successfully, the name in db_name, model_db_path and covariance must be consistent. Which means that if we want to have a trait called 50_raw to perform the analysis, we must have corresponding db file 50_raw.db in model_db_path folder and covariance file 50_raw_cov.txt.gz and 50_raw_cov_tn.txt.gz in covariance folder.

## Multi-variate analysis

In order to eliminate effect between traits, users can also perform multi-variate analysis with more than one traits as input, the command will be similar to the command above. The only difference is instead of using python BADGERS.py in the first line, user needs to use python BADGERS_mult.py

By this, all traits in db_name file will be involved in multi-variate analysis. We encourage to first perform selection process like clustering and make trait number involved in multi-variate analysis to be less than 100 for better statistical power. (Since more than 100 traits will make imputed inverse covariance matrix to be inaccurate)

## Acknowledgement
Part of the code is modified from MetaXcan https://github.com/hakyimlab/MetaXcan. We thank the authors for sharing the code.

## Reference
**Yan et al. (2018). Biobank-wide association scan identifies risk factors for late-onset Alzheimer’s disease and endophenotypes. bioRxiv, 468306.** [Link](https://www.biorxiv.org/content/early/2018/11/12/468306)
