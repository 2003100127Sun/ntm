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

* initial data for a backcross
```
python bcNTMapping -fg geno.csv -fp pheno.csv
```
> parameters illustration
```
-fg --filegeno -> a file for phenotype. required.
-fp --filepheno -> a file for genotype which can consist of different chromosomes. required
```

```
python bcNTMapping -fg geno.csv -fp pheno.csv -spt 1
```
> parameters illustration
```
-spt --split -> split your genotype file into different pieces of chromosome group each of which has the information about entire intervals, Morgan distances and markers.
-nt -- NT values using numerical trajectory measurement
```
> file format
1. phenotype file: The file contains growth data observed over time, headers of the data, with the columns for time and rows for individuals.
2. genotype file: The file contains genomic distances, chromosome groups and markers. Note that the markers should be . In a backcross, the format of the phenotype file should 0 for aa and 1 for Aa and -1 for missing data. Or, you can specify the genotype of flanking mark with aa and Aa. If the number form (0,1,-1) is used, the final form used for QTL mapping will be performed by a built-in method 'split'.

* obtain NT values

```
python bcNTMapping -fg geno.csv -fp pheno.csv -nt 1
```
> parameters illustration
```
-nt -- NT values using numerical trajectory measurement
```
> explanation
return a 1d array. Each element represents a growth trajectory for individual i changing over time.