Crown Tutorial
================

## 1. Overview
This repository provides the code needed to compute the CROWN method for **C**onstrained **R**esidual Nodewise **O**ptimal **W**eight Regressio**n**, which can be used to construct constrained portfolios in a high-dimensional scenario (p>T). It also shows different methods in obtaining a solution for high dimensional portfolio with a large number of assets. Specifically, we consider portfolios with tracking error constraints, portfolios with tracking error jointly with weight (equality or inequality) restrictions, and portfolios with only weight restrictions. 

Firstly, we need an estimator of precision matrix for the return data.

This repository shows the method **CROWN** along with 4 other popular methods:
1. **NODEWISE** from [A Nodewise Regression Approach to Estimating Large Portfolios]
2. **POET** from [Large Covariance Estimation by Thresholding Principal Orthogonal Complements]
3. **NLS**  from [Nonlinear Shrinkage of the Covariance Matrix for Portfolio Selection: Markowitz Meets Goldilocks]
4. **SF-NLS** from [Nonlinear Shrinkage of the Covariance Matrix for Portfolio Selection: Markowitz Meets Goldilocks]

The paper detailing the theory and methodology behind the CROWN estimator is publicly available and can be accessed at [Navigating Complexity: Constrained Portfolio Analysis in High Dimensions with Tracking Error and Weight Constraints](https://arxiv.org/abs/2402.17523).


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


| Method | TE      | Weight_ER | Risk_ER | SR_ER  |
|--------|---------|-----------|---------|--------|
| NW     | 1.2040  | 8.7975    | 41.4364 | **0.9200** |
| CROWN  | **0.1318**  | **1.5408**    | **0.1008**  | 0.9996 |
| POET   | 0.3382  | 5.0783    | 6.6797  | 0.9917 |
| NLS    | 0.1349  | 4.3591    | 4.1078  | 0.9981 |
| SFNL   | 0.1452  | 2.4806    | 0.7936  | 0.9987 |


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

## References

Callot, L., M. Caner, O. Onder, and E. Ulasan (2021). A nodewise regression approach to estimating large portfolios. Journal of Business and Economic Statistics 39, 520–531.

Caner, M., Fan, Q., & Li, Y. (2024). Navigating Complexity: Constrained Portfolio Analysis in High Dimensions with Tracking Error and Weight Constraints. arXiv preprint arXiv:2402.17523.

Fan, J., Y. Liao, and M. Mincheva (2013). Large covariance estimation by thresholding principal orthogonal complements. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 75 (4), 603–680.

Ledoit, O, M. and M. Wolf (2017). Nonlinear shrinkage of the covariance matrix for portfolio selection: Markowitz meets goldilocks. Review of Financial Studies 30, 4349–4388.
