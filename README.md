# BADGERS

BADGERS(Biobank-wide Association Discovery using GEnetic Risk Scores), a powerful method to perform polygenic score-based biobank-wide association scan.

## Prerequisites

The software is developed and tested in Linux and Mac OS environments.
- Python 2.7
- numpy
- scipy
- pandas

## Quick Start 

This section demonstrates the fundamental usage of BADGERS, which is to conduct the association between selelcted UK-biobank traits (total of 4357) and a complex disease sumstats. If you want to run BADGERS with other complex traits. You need to create corresponding [weight database](https://github.com/qlu-lab/BADGERS/wiki/Create-db-files) and [covariance](https://github.com/qlu-lab/BADGERS/wiki/Create-covariance-file) to replace UK-biobank traits as input in the sample below.

### Step1: Downloads BADGERS

```
$ git clone https://github.com/qlu-lab/BADGERS
$ cd ./BADGERS
```
In order to test whether BADGERS function correctly in your machine you can run [sample test](https://github.com/qlu-lab/BADGERS/wiki/Sample-test) for it

### Step2: Downloads UK-biobank traits weight database and covariance
```
$ wget ftp://ftp.biostat.wisc.edu/pub/lu_group/BADGERS/UK_biobank_input
$ tar xf UK_biobank_input.tar.gz
```
This folder will include the following files/folders:
```
cov/      ## covariance file for traits
weight_db/      ## weight database for traits
UKbiobank_4357_inputlist.csv      ## name for traits
```
### Step3: Perform analysis

Here is how you get the association between UK_biobank traits and complex disease

```
python GAACT.py \
--model_db_path /weight_db \
--covariance /cov \
--gwas_path IGAP.txt \
--snp_column MarkerName \
--effect_allele_column Allele1 \
--non_effect_allele_column Allele2 \
--pvalue_column P.value \
--beta_column Effect \
--output_file output.csv \
--db_name selected_traits.csv 
```
where
- *--model_db_path:*

    The location of database files for all traits.

- *--covariance:*

    The location of covariance files for all traits.

- *--gwas_path:*

    The location of input GWAS data(user input).

- *--snp_column:*

    Name of column holding SNP data.

- *--effect_allele_column:*

    Name of column holding effect allele data.

- *--non_effect_allele_column:*

    Name of column holding "other/non effect" allele data.

- *--pvalue_column:*

    Name of column holding p-values data.

- *--beta_column:*

    Name of column holding beta data.

- *output_file:*

    Location where results will be saved.

- *db_name:*

    Database file name for traits you want to perform analysis, for example if you want to test the association between 50_raw.db,             20016_raw.db and input disease, for db_name input you should have a selected_traits.csv file like:
<pre>
    >   50_raw.db,    Standing height\
    >   20016_raw.db, Fluid intelligence score
<pre\>
## Acknowledgement
Part of the code is modified from MetaXcan https://github.com/hakyimlab/MetaXcan. We thank the authors for sharing the code.

## Reference
