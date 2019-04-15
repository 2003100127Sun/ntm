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
	* Here is an option for you. You use R if installing R package qtl to allow for a fast speed and better performance [recommend option] or you can use a self-programming module of mixed model of interval mapping for a BC otherwise. 

## Usage

1. initial data for a backcross
	```
	python bcNTMapping -fg geno.csv -fp pheno.csv
	```
	> parameters illustration
	```
	-fg --filegeno -> a file for phenotype. required.
	-fp --filepheno -> a file for genotype which can consist of different chromosomes. required
	```

2. split and construct data
	```
	python bcNTMapping -fg geno.csv -fp pheno.csv -spt 1
	```
	> parameters illustration
	```
	-spt --split -> split your genotype file into different pieces of chromosome group each of which has the information about entire intervals, Morgan distances and markers.
	```
	> file format
	1. phenotype file: The file contains growth data observed over time, headers of the data, with the columns for time and rows for individuals.
	2. genotype file: The file contains genomic distances, chromosome groups and markers. Note that the markers should be . In a backcross, the format of the phenotype file should 0 for aa and 1 for Aa and -1 for missing data. Or, you can specify the genotype of flanking mark with aa and Aa. If the number form (0,1,-1) is used, the final form used for QTL mapping will be performed by a built-in method 'split'.

3. obtain NT values
	NT values are calculated with the divided difference.
	```
	python bcNTMapping -fg geno.csv -fp pheno.csv -nt 1
	```
	> parameters illustration
	```
	-nt -- NT values -> using numerical trajectory measurement ***
	```
	> explanation: 
	return a 1d array. Each element represents a growth trajectory for individual i changing over time.

4. parser data
	* profile of interval table 
	```
	python bcNTMapping -fg geno.csv -fp pheno.csv -it 1
	```
	> parameters illustration
	```
	-it -- itable -> interval table
	```
	> explanation: 
	return a 1d array. Each element represents a growth trajectory for individual i changing over time.
	* profile of intervals, genomic regions, groups of chromosome.
	```
	python bcNTMapping -fg geno.csv -fp pheno.csv -bpf 1
	```
	> parameters illustration
	```
	-bpf --basicprofile -> basic profile of genotypic info
	```
	> explanation: 
	return a 1d array. Each element represents intervals, genomic regions, groups of chromosome.

5. permutation test
	```python batchPermutation -ng 22 - np 1000```
	> parameters illustration
	```
	-ng --ngroups -> number of chromosome groups. default value, 1.
	-np --npermutation -> number of permutation test. default value, 1000.
	```