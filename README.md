# NT_Measurement - Backbone Release
![](https://img.shields.io/badge/NT_Measurement-Backbone-519dd9.svg)
![](https://img.shields.io/badge/last_released_date-April_2019-green.svg)
## Overview
##### The repository will be updated later for the backbone part of NT_Measurement.

* stamp 1
    + version 1.0
    + time: April 14, 2019

## Installation
* install necessary package
```
pip install -r requirements.txt
```

## Usage
* permutation test
```python batchPermutation -ng 22 - np 1000```
> parameters illustration
```
-ng --help number of chromosome groups. default value, 1.
-np --help number of permutation test. default value, 1000.
```

* parser data for a backcross
```
python bcNTMapping -fg geno.csv -fp pheno.csv
```
> parameters illustration
```
-fg --help a file for phenotype. required.
-fp --help a file for genotype which can consist of different chromosomes. required
```
> file format
	1. phenotype file: a file contains growth data observed over time, headers of the data, with the columns for time and rows for individuals.
	2. genotype file: a file contains genomic distance, chromosome group and markers.

