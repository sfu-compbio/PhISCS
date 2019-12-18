# PhISCS

PhISCS is a tool for sub-perfect tumor phylogeny reconstruction via integrative use of single-cell and bulk sequencing data. If bulk sequencing data is used, we expect that mutations originate from diploid regions of the genome. Due to variance in VAF values, we recommend the use of bulk data in cases when sequencing depth is at least 1000x (haploid coverage). As output, PhISCS reports tree of tumor evolution together with a set of eliminated mutations, where eliminated mutations represent mutations violating Infinite Sites Assumption (due to deletion of variant allele or due to recurrent mutation) or mutations affected by copy number aberrations that were missed during the tumor copy number profiling (e.g. gain of non-variant allele).

PhISCS has been published in Genome Research [(doi:10.1101/gr.234435.118)](http://www.genome.org/cgi/doi/10.1101/gr.234435.118)

## Contents
  1. [Installation](#installation)
     * [PhISCS-I](#installationilp)
       * [Prerequisite: ILP solver](#prerequisiteilp)
     * [PhISCS-B](#installationcsp)
       * [Prerequisite: CSP solver](#prerequisitecsp)
  2. [Running](#running)
     * [Input](#input)
       * [Single-cell Matrix](#singlecellmatrix)
       * [Bulk Data](#bulkdata)
     * [Output](#output)
       * [Log File](#logfile)
       * [Output Matrix File](#outputmatrixfile)
     * [Parameters](#parameters)
  3. [Example](#example)
  4. [Contact](#contact)

<a name="installation"></a>
## Installation
PhISCS is written in Python and C. It supports both Python 2.7 and 3. Currently it is intended to be run on POSIX-based systems (only Linux and macOS have been tested).

> **RECOMENDATION**:  At the moment, in cases when both, single-cell and bulk data are used as input, we recommend the use of PhISCS-I over PhISCS-B (due to more thorough tests and software validation that we have performed for PhISCS-I). However, when single-cell data is the only input, we have extensively tested both implementations and, since PhISCS-B can have potential running time advantage in this case, we recommend its use over PhISCS-I.   

<a name="installationilp"></a>
### PhISCS-I
```
git clone --recursive https://github.com/sfu-compbio/PhISCS.git
cd PhISCS
python PhISCS-I --help
```

<a name="prerequisiteilp"></a>
#### Prerequisite: ILP solver

In order to run PhISCS-I, the main requirement is the installation of Gurobi solver. [Gurobi](http://www.gurobi.com) a commercial solver which is free for academic purposes. After installing it, installation of `gurobipy` package is necessary prior to being able to successfully run PhISCS-I (below we provide some examples of the input and commands used to run the tool).


<a name="installationcsp"></a>
### PhISCS-B

```
git clone --recursive https://github.com/sfu-compbio/PhISCS.git
cd PhISCS
./PhISCS-B-configure
python PhISCS-B --help
```

<a name="prerequisitecsp"></a>
#### Prerequisite: CSP solver

Some of CSP solver have been already included in the PhISCS package. There is an option to add a new CSP solver to PhISCS-B by provinding a path to the exe file of the desired CSP solver.


<a name="running"></a>
## Running

<a name="input"></a>
### Input

<a name="singlecellmatrix"></a>
#### 1. Single-cell Matrix
Single-cell input is assumed to be represented in the form of ternary, __tab-delimited__, matrix with rows corresponding to single-cells and columns corresponding to mutations. We assume that this file contains headers and that matrix is ternary matrix with 0 denoting the absence and 1 denoting the presence of mutation in a given cell, whereas ? represents the lack of information about presence/absence of mutation in a given cell (i.e. missing entry). __In order to simplify parsing of the matrix, we also assume that upper left corner equals to string `cellID/mutID`__.

Below is an example of single-cell data matrix. Note that mutation and cell names are arbitrary strings not containing tabs or spaces, however they must be unique.
```
cellID/mutID  mut0  mut1  mut2  mut3  mut4  mut5  mut6  mut7
cell0         0     0     ?     0     0     0     0     0
cell1         0     ?     1     0     0     0     1     1
cell2         0     0     1     0     0     0     1     1
cell3         1     1     0     0     0     0     0     0
cell4         0     0     1     0     0     0     0     0
cell5         1     0     0     0     0     0     0     0
cell6         0     0     1     0     0     0     1     1
cell7         0     0     1     0     0     0     0     0
cell8         ?     0     0     0     ?     0     ?     1
cell9         0     1     0     0     0     0     0     0
```

<a name="bulkdata"></a>
#### 2. Bulk Data
As bulk data input, we also expect __tab-delimited__ file with the following columns:
  
**ID** which represents mutational ID (used in single-cell data matrix for the same mutation)  
**Chromosome** which represents chromosome of the mutation (any string not containing tabs or empty spaces)  
**Position** which represents position (on chromosome) of the mutation (any string/number not containing tabs or empty spaces)  
**MutantCount** is the number of mutant reads in the bulk data. If multiple bulk samples are used, values are semicolon-delimited and provided in the sorted order of samples (this order is expected to be same for all mutations, e.g. first number always representing read count in sample 1, second number in sample 2 etc.)  
**ReferenceCount** is the number of reference reads in the bulk data. If multiple bulk samples are used, values are semicolon-delimited and provided in the sorted order of samples (this order is expected to be same for all mutations, e.g. first number always representing read count in sample 1, second number in sample 2 etc.)  
**INFO** which contains additional information about the mutation and is semicolon-delimited. Entries in this column are of the form: entryID=values, where values are delimited by commas. An example of INFO column is:
"sampleIDs=S0,S1;synonymous=false;exonic=true". The only obligatory information required now is information about sample origins (in cases of absence of them, arbitrary distinct strings can be used, e.g. sampleIDs=S0,S1,S2;)



As an example:
```
ID    Chromosome  Position  MutantCount     ReferenceCount    INFO
mut0  1           0         766;511;688     4234;4489;4312    sampleIDs=primary,metastasis1,metastasis2
mut1  1           1         719;479;719     4281;4521;4281    sampleIDs=primary,metastasis1,metastasis2
mut2  1           2         1246;1094;859   3754;3906;4141    sampleIDs=primary,metastasis1,metastasis2
mut3  1           3         298;226;272     4702;4774;4728    sampleIDs=primary,metastasis1,metastasis2
mut4  1           4         353;227;255     4647;4773;4745    sampleIDs=primary,metastasis1,metastasis2
mut5  1           5         306;232;279     4694;4768;4721    sampleIDs=primary,metastasis1,metastasis2
mut6  1           6         725;449;492     4275;4551;4508    sampleIDs=primary,metastasis1,metastasis2
mut7  1           7         703;417;507     4297;4583;4493    sampleIDs=primary,metastasis1,metastasis2

```

(in the example of bulk file shown above, we have that for mut0 number of mutant and reference reads in the first sample are respectively 766 and 4234, in the second sample 511 and 4489 and in the third sample 688 and 4312).

<a name="output"></a>
### Output
The program will generate two files in **OUT_DIR** folder (which is set by argument -o or --outDir). This folder will be created automatically if it does not exist.

<a name="outputmatrixfile"></a>
#### 1. Output Matrix File
The output matrix is also a tab-delimited file having the same format as the input matrix, except that eliminated mutations (columns) are excluded (so, in case when mutation elimination is allowed, this matrix typically contains less columns than the input matrix). Output matrix represents genotypes-corrected matrix (where false positives and false negatives from the input are corrected and each of the missing entries set to 0 or 1). Suppose the input file is **INPUT_MATRIX.ext**, the output matrix will be stored in file **OUT_DIR/INPUT_MATRIX.CFMatrix**. For example:
```
 input file: data/ALL2.SC
output file: OUT_DIR/ALL2.CFMatrix
```

<a name="logfile"></a>
#### 2. Log File
Log file contains various information about the particular run of PhISCS (e.g. eliminated mutations or likelihood value). The interpretation of the relevant reported entries in this file is self-evident. Suppose the input file is **INPUT_MATRIX.ext**, the log will be stored in file **OUT_DIR/INPUT_MATRIX.log**. For example:
```
input file: data/ALL2.SC
  log file: OUT_DIR/ALL2.log
```

<a name="parameters"></a>
### Parameters
| Parameter  | Description                                                                                | Default  | Mandatory      |
|------------|--------------------------------------------------------------------------------------------|----------|----------------|
| -SCFile    | Path to single-cell data matrix file                                                       | -        | :radio_button: |
| -fn        | Probablity of false negative                                                               | -        | :radio_button: |
| -fp        | Probablity of false positive                                                               | -        | :radio_button: |
| -o         | Output directory                                                                           | current  | :white_circle: |
| -kmax      | Max number of mutations to be eliminated                                                   | 0        | :white_circle: |
| -threads   | Number of threads (supported by PhISCS-I)                                                  | 1        | :white_circle: |
| -bulkFile  | Path to bulk data file                                                                     | -        | :white_circle: |
| -delta     | Delta parameter accounting for VAF variance                                                | 0.20     | :white_circle: |
| -time      | Max time (in seconds) allowed for the computation                                          | 24 hours | :white_circle: |
| --drawTree | Draw output tree with Graphviz                                                             | -        | :white_circle: |

<a name="example"></a>
## Example

For running PhISCS without VAFs information and without ISA violations:
```
python PhISCS-I -SCFile example/input.SC -fn 0.2 -fp 0.0001 -o result/
```

For running PhISCS without VAFs information but with ISA violations:
```
python PhISCS-I -SCFile example/input.SC -fn 0.2 -fp 0.0001 -o result/ -kmax 1
```

For running PhISCS with both VAFs information and ISA violations (with time limit of 24 hours):
```
python PhISCS-I -SCFile example/input.SC -fn 0.2 -fp 0.0001 -o result/ -kmax 1 -bulkFile example/input.bulk -time 86400
```

For running PhISCS with VAFs information but no ISA violations (with drawing the output tree):
```
python PhISCS-I -SCFile example/input.SC -fn 0.2 -fp 0.0001 -o result/ -bulkFile example/input.bulk --drawTree
```

<a name="contact"></a>
## Contact
If you have any questions please e-mail us at smalikic@sfu.ca or frashidi@iu.edu.
