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

### Step2: Downloads 
