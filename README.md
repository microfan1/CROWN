Crown Tutorial
================

## 1. Overview
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


## 2. Required packages
You need to download matlab toolbox RunRCode for R functions calling.

You can install the required package with

```
install.pacakges(c("R.matlab, glmnet, nlshrink, PDSCE, iterators"))
```

## 3. Path Editing

Please modify the working path in the R language function, the following is the relevant file path:

```
\code\methods\NODEWISE\Nodewise_Run_R.R

\code\methods\POET\POET_R.R

\code\methods\SF-NLS\SF-NLS_Run_R.R
```

## 4. Repository Contents
This repository contains the Matlab code and R scripts used for the comparison of 5 methods. 

- Methods/
   - CROWN: 
     - glmnet_matlab: Package used by `crown.m`.
     - `crown.m`: Function to get crown estimator.
   - NLS:
     - `nls_covMarket.m`: Function to get nls estimator.
   - NODEWISE: 
     - `Nodewise_Run_R.R`: Function to get Nodewise estimator.
     - `Nodewise_source.R`: Functions used by `Nodewise_Run_R.R`.
   - POET:
     - POET: Package used by `POET_R.R`.
     - `POET.m`: Code to link R and Matlab.
     - `POET_R.R`: Code to get POET estimator.
   - SF-NLS:
     - `SF-NLS_Run_R.R`: Code to get SF-NLS estimator.
## 5. Example of Usage
`Example.mlx` shows how to use crown and 4 baseline methods for portfolios with 4 kinds of constraints. They are portfolios with tracking error constraints, portfolios with tracking error jointly with weight (equality or inequality) restrictions, and portfolios with only weight restrictions. 

Please modify the working path and Rpath in `Example.mlx`, `Nodewise_Run_R.R`, `POET_R.R` and `SF-NLS_Run_R.R` first.

### 5.1 Example and Data Loading: 
1. 180 stocks with 150 trading records;
2. Let Tracing Error constraint equals 0.2;
3. Let sum of first 10 weights equals to 0.2 for weights constraint.
        
```
% Set your portfolio and sample size
pNpairs = [180,150];p=pNpairs(1);N=pNpairs(2);

% Set your Tracing Error Constraint
TEC = 0.2;

% Set your Weights Constraint
oneR = [ones(10,1);zeros(p-10,1)];
omega = 0.2;
```
Change parameters to suit your demand. For data loading below, feel free to change them if you want:

```
%3-Factor
flag=3;k = 3;
%Data Loading 
beta_total=mvnrnd([0.001,0.001,0.001],0.1.*eye(3),p);
bf = [0.005,0,0;0,0.005,0;0,0,0.005];% B Ft + et
mu_f=[0.005,0.005,0.005]';%kx1
cov_f = eye(3);
```
Use `Data_Generation.m` to generate simulation data based on these loading. 

### 5.2 Estimation of 5 methods
Use these code to estimate precision matrix:
```
%NODEWISE: 
  save(fullfile(pwd, '\methods\NODEWISE\to_Nodewise_R.mat'),'R','fac')
  RunRcode(fullfile(pwd,'\methods\NODEWISE\Nodewise_Run_R.R'),Rpath);
  inv_Sigma_nw = csvread(fullfile(pwd,'\methods\NODEWISE\inv_Sigma_nw.csv'));

%CROWN
  addpath methods\CROWN
  addpath methods\CROWN\glmnet_matlab
  O=crown(R);
  inv_Sigma_crown=O-O*beta/(eye(k)+beta'*(O+O')/2*beta)*beta'*O;

%POET
  addpath methods\POET
  Sigma_poet = POET(R,fullfile(pwd,'\methods\POET'),Rpath);
  inv_Sigma_poet=inv(Sigma_poet);

%NLS
  addpath methods\NLS
  Sigma_nls=nls_covMarket(R);
  inv_Sigma_nls=inv(Sigma_nls);

%SF-NLS
  save(fullfile(pwd, '\methods\SF-NLS\to_SF-NLS_R.mat'),'R','fac')
  RunRcode(fullfile(pwd,'\methods\SF-NLS\SF-NLS_Run_R.R'),Rpath);
  Sigma_sfnls = csvread(fullfile(pwd,'\methods\SF-NLS\Sigma_SF-NLS.csv'));
  inv_Sigma_sfnls = inv(Sigma_sfnls);
```

### 5.3 Get Weights and Performance
1. Only Consider Tracing Error Constraint
   
| Method | TE     | Weight_ER | Risk_ER   | SR_ER   |
|--------|--------|-----------|-----------|---------|
| NW     | 2.6408 | 6.5044    | 182.9501  | 0.0103  |
| CROWN  | **0.2284** | **0.1808**    | **0.0431**    | **0.0082**  |
| POET   | 0.5542 | 4.2305    | 16.1620   | 0.3471  |
| NLS    | 0.3482 | 3.8828    | 12.4516   | 0.2734  |
| SF-NLS | 0.3596 | 2.4863    | 7.1400    | 0.1972  |

2. Only consider Weights Constraint

| Method | TE     | Weight_ER | Risk_ER | SR_ER  |
|--------|--------|-----------|---------|--------|
| NW     | **0.0726** | 1.1314    | 0.5368  | **0.9272** |
| CROWN  | 0.0738 | **1.0465**    | 0.6625  | 0.9996 |
| POET   | 0.0774 | 1.2298    | 0.5734  | 0.9920 |
| NLS    | 0.0881 | 1.3566    | **0.5012**  | 0.9985 |
| SFNL   | 0.1024 | 1.1455    | 0.5931  | 0.9989 |


3. Tracing Error + Weights Constraint

| Method | TE     | Weight_ER | Risk_ER   | SR_ER   |
|--------|--------|-----------|-----------|---------|
| NW     | 2.6124 | 6.5295    | 178.9517  | **0.0120**  |
| CROWN  | **0.2431** | **0.4723**    | **0.1352**    | 0.1752  |
| POET   | 0.5516 | 4.2925    | 15.9899   | 0.3491  |
| NLS    | 0.3493 | 3.9574    | 12.7060   | 0.2872  |
| SFNL   | 0.3677 | 2.5603    | 7.3275    | 0.2161  |


4. Tracing Error +  Weights Constraint though non-binding 

| Method | TE     | Weight_ER | Risk_ER   | SR_ER   |
|--------|--------|-----------|-----------|---------|
| NW     | 2.6408 | 6.5044    | 182.9501  | 0.0103  |
| CROWN  | **0.2284** | **0.1808**    | **0.0431**    | **0.0082**  |
| POET   | 0.5542 | 4.2305    | 16.1620   | 0.3471  |
| NLS    | 0.3482 | 3.8828    | 12.4516   | 0.2734  |
| SFNL   | 0.3596 | 2.4863    | 7.1400    | 0.1972  |


