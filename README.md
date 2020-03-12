# Co-expression analysis in Go
An implementation of weighted gene co-expression network analysis in Go. This is a final project for CMU's Programming for Scientists class. 

#### Authors
* **Serena Abraham** - *matrices.go*
* **Jon Luo** - *clustering.go*
* **Hanxi Xiao** - *io.go*

All authors contributed to each set.

Example data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115522

## Getting Started

### Prerequisites
You will need ~10GB free RAM.

Gonum is required to run the project. 
```
go get -u gonum.org/v1/gonum/...
```
Unzip data.zip before running.

### Install and run

```
go build
unzip data.zip
./coexpression
```

You will be prompted to choose testing options. We recommend Pearsons and Signed correlation network.

----------------------------------------
CO-EXPRESSION GRAPHS AND THEIR ANALYSIS
----------------------------------------
For the covariance tests, please enter P for Pearsons (faster) or B for BiWeightedCorrelation:

Enter S to build a Signed correlation network or U for Unsigned:

Gene clusters will be saved to clusters.txt with one cluster per line.

Approximate runtime on a Ryzen 3700X 8-core/16-thread CPU is 15 minutes.
