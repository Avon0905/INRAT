# INRAT
A computational method for inferring gene regulatory network by integrating time series single-cell RNA-seq and ATAC-seq data. Taking single cell expression data and epigenetic data as input, INRAT identifies regulatory relationships between target genes and transcription factors by using non-linear ODE model. 

## Requirements
INRAT is written with R.

## Download
```
git clone https://github.com/Avon0905/INRAT
cd INRAT
```

## Implementation
```
Rscript run_INRAT.R <Input_file1> <Input_file2> <Input_file3> <Output_dir> <G> <D> <C> <I>
```
* Input_file1 : G x C matrix of expression data
* Input_file2 : Time point data (e.g. pseudo-time data)
* Input_file3 : Binding Score file from ATAC-seq data
* Output_dir : Result files are outputted in this directory
* G : The number of transcription factors
* D : The number of z
* C : The number of cells
* I : The number of iterations of optimization


##### Format of Input_file1
The Input_file1 is the G x C matrix of expression data (separated with 'TAB').
Each row corresponds to each gene, and each column corresponds to each cell.

##### Example of Input_file1
```
1.24	1.21	1.28	...
0.0 	0.19	0.0	...
.
.
.
```

##### Format of Input_file2
The Input_file2 contains the time point data (pseudo-time) of each cell.

* Col1 : Information of a cell (e.g. index of a cell, experimental time point)
* Col2 : Time parameter (e.g. pseudo-time) (normalized from 0.0 to 1.0)

##### Example of Input_file2
```
0	0.065
0	0.037
0	0.007
.
.
.
72	0.873
72	0.964
```

