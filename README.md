Crown Tutorial
================

## Overview
This repository shows comparison of different methods in obtaining a solution for high dimensional portfolio with a large number of assets. That is, we consider portfolios with tracking error constraints, portfolios with tracking error jointly with weight (equality or inequality) restrictions, and portfolios with only weight restrictions. 

The key of this process is to get a good estimator of precision matrix for the return data.

This repository compares our method **CROWN** with 4 baseline methods:
1. **NODEWISE** from [A Nodewise Regression Approach to Estimating Large Portfolios](https://arxiv.org/pdf/1611.07347)
2. **POET** from [Large Covariance Estimation by Thresholding Principal Orthogonal Complements](https://arxiv.org/pdf/1201.0175)
3. **NLS**  from [Nonlinear Shrinkage of the Covariance
Matrix for Portfolio Selection: Markowitz
Meets Goldilocks](http://www.ledoit.net/Goldilocks_RFS_2017.pdf)
4. **SF-NLS** from [Nonlinear Shrinkage of the Covariance
Matrix for Portfolio Selection: Markowitz
Meets Goldilocks](http://www.ledoit.net/Goldilocks_RFS_2017.pdf)


## Required packages
You need to download matlab toolbox RunRCode for R functions calling.

You can install the required package with

```
install.pacakges(c("R.matlab, glmnet, nlshrink, PDSCE, iterators"))
```

## Path Editing

Please modify the working path and Rpath in the matlab file:

```
pwd = 'E:\elective\RA\CUHK\Tutorial\code';cd(pwd);

Rpath = 'E:\major\IT\R\R-4.3.3\bin';
```

Please modify the working path in the R language function, the following is the relevant file path:

```
\code\methods\NODEWISE\Nodewise_Run_R.R

\code\methods\POET\POET_R.R

\code\methods\SF-NLS\SF-NLS_Run_R.R
```

## Repository Contents

## Example of Usage

## References

