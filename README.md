Crown
================

## 1. Overview
This repository provides the code needed to compute the CROWN method for **C**onstrained **R**esidual Nodewise **O**ptimal **W**eight Regressio**n**, which can be used to construct constrained portfolios in a high-dimensional scenario (p>T). It also shows different methods in obtaining a solution for high dimensional portfolio with a large number of assets. Specifically, we consider portfolios with tracking error constraints, portfolios with tracking error jointly with weight (equality or inequality) restrictions, and portfolios with only weight restrictions. 

Firstly, we need an estimator of precision matrix for the return data.

This repository shows the method **CROWN** using a residual-based nodewise regression to get the estimate of covariance, along with 4 other popular methods:
1. **NODEWISE** from [A Nodewise Regression Approach to Estimating Large Portfolios]
2. **POET** from [Large Covariance Estimation by Thresholding Principal Orthogonal Complements]
3. **NLS**  from [Nonlinear Shrinkage of the Covariance Matrix for Portfolio Selection: Markowitz Meets Goldilocks]
4. **SF-NLS** from [Nonlinear Shrinkage of the Covariance Matrix for Portfolio Selection: Markowitz Meets Goldilocks]

The paper detailing the methodology for the CROWN estimator is available at [Navigating Complexity: Constrained Portfolio Analysis in High Dimensions with Tracking Error and Weight Constraints](https://arxiv.org/abs/2402.17523).


## 2. Required packages
You need to download matlab toolbox RunRCode for R functions calling.

You can install the required package with

```
install.pacakges(c("R.matlab, glmnet, nlshrink, PDSCE, iterators"))
```

## 3. Path Editing

Please modify the working path in the R language function, the following is the relevant file path:

```
\code\covariance estimation\NODEWISE\Nodewise_Run_R.R

\code\covariance estimation\POET\POET_R.R

\code\covariance estimation\SF-NLS\SF-NLS_Run_R.R
```

## 4. Repository Contents
This repository contains the Matlab code and R scripts used for the comparison of 5 methods. 

- covariance estimation/
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
```
For data loading below, feel free to change them (or load your own data in similar format):
```
### 5.1 Example and Data Loading: 
1. 180 stocks with 150 trading records;
2. 3 factors with 150 observations;
3. Let Tracking Error constraint equals 0.2;
4. Let sum of first 10 weights equals to 0.2 for weights constraint.
        
```
% Set your portfolio and sample size
pNpairs = [180,150];p=pNpairs(1);N=pNpairs(2);

% Set your Tracking Error Constraint
TEC = 0.2;

% Set your Weights Constraint
oneR = [ones(10,1);zeros(p-10,1)];
omega = 0.2;
```
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
```
### 5.2 Estimation of 5 methods
Use these code to estimate precision matrix:
```
%NODEWISE: 
  save(fullfile(pwd, '\covariance estimation\NODEWISE\to_Nodewise_R.mat'),'R','fac')
  RunRcode(fullfile(pwd,'\covariance estimation\NODEWISE\Nodewise_Run_R.R'),Rpath);
  inv_Sigma_nw = csvread(fullfile(pwd,'\covariance estimation\NODEWISE\inv_Sigma_nw.csv'));

%CROWN
  addpath('covariance estimation\CROWN')
  addpath('covariance estimation\CROWN\glmnet_matlab')
  O=crown(R);
  inv_Sigma_crown=O-O*beta/(eye(k)+beta'*(O+O')/2*beta)*beta'*O;

%POET
  addpath('covariance estimation\POET')
  Sigma_poet = POET(R,fullfile(pwd,'\covariance estimation\POET'),Rpath);
  inv_Sigma_poet=inv(Sigma_poet);

%NLS
  addpath('covariance estimation\NLS')
  Sigma_nls=nls_covMarket(R);
  inv_Sigma_nls=inv(Sigma_nls);

%SF-NLS
  save(fullfile(pwd, '\covariance estimation\SF-NLS\to_SF-NLS_R.mat'),'R','fac')
  RunRcode(fullfile(pwd,'\covariance estimation\SF-NLS\SF-NLS_Run_R.R'),Rpath);
  Sigma_sfnls = csvread(fullfile(pwd,'\covariance estimation\SF-NLS\Sigma_SF-NLS.csv'));
  inv_Sigma_sfnls = inv(Sigma_sfnls);
```
```
### 5.3 Get Weights and Performance
It's worth mentioning that the TE constrained weights we get fron 'Get_weights' functions are the differences between the studied portfolio and a benchmark index, which means they should be added by tracking index b to get the real weights.
1. Only Consider Tracking Error Constraint
   
| Method | TE     | Weight_ER | Risk_ER   | SR_ER   |
|--------|--------|-----------|-----------|---------|
| NW     | 39.3776| 97.5653   | 182.7355  | 0.0103  |
| CROWN  | **2.9644** | **2.3114**    | **0.0128**    | **0.0076**  |
| POET   | 7.7450 | 63.4571   | 16.0826   | 0.3471  |
| NLS    | 4.8113 | 58.2418   | 12.3971   | 0.2741  |
| SF-NLS | 4.9803 | 37.2946   | 7.0945    | 0.1976  |

2. Only consider Weights Constraint


| Method | TE        | Weight_ER | Risk_ER  | SR_ER    |
|--------|-----------|-----------|-----------|----------|
| NW     | 36.8717   | 98.0123   | 36237.00  | **0.0103**   |
| CROWN  | **2.9948**    | **11.5354**   | **233.6463**  | 0.0771   |
| POET   | 7.2979    | 61.7397   | 3403.20   | 0.3470   |
| NLS    | 4.5391    | 56.9357   | 2678.30   | 0.2741   |
| SFNL   | 4.7034    | 37.5961   | 1620.80   | 0.1983   |


3. Tracking Error + Weights Constraint

| Method | TE     | Weight_ER | Risk_ER   | SR_ER   |
|--------|--------|-----------|-----------|---------|
| NW     | 39.2824| 97.3143   | 182.0070  | 0.0106  |
| CROWN  | **2.9586** | **2.4644**    | **0.0124**    | **0.0087**  |
| POET   | 7.7380 | 63.2637   | 16.0354   | 0.3458  |
| NLS    | 4.8114 | 58.1513   | 12.3885   | 0.2738  |
| SFNL   | 4.9807 | 37.2932   | 7.0947    | 0.1976  |


4. Tracking Error +  Weights Constraint though non-binding 

| Method | TE     | Weight_ER | Risk_ER   | SR_ER   |
|--------|--------|-----------|-----------|---------|
| NW     | 39.3776| 97.5653   | 182.7355  | 0.0103  |
| CROWN  | **2.9644** | **2.3114**    | **0.0128**    | **0.0076**  |
| POET   | 7.7450 | 63.4571   | 16.0826   | 0.3471  |
| NLS    | 4.8113 | 58.2418   | 12.3971   | 0.2741  |
| SF-NLS | 4.9803 | 37.2946   | 7.0945    | 0.1976  |

## References

Callot, L., M. Caner, O. Onder, and E. Ulasan (2021). A nodewise regression approach to estimating large portfolios. Journal of Business and Economic Statistics 39, 520–531.

Caner, M., Fan, Q., and Li, Y. (2024). Navigating Complexity: Constrained Portfolio Analysis in High Dimensions with Tracking Error and Weight Constraints. arXiv preprint arXiv:2402.17523.

Caner, M., Medeiros, M., and G. Vasconcelos (2023). Sharpe Ratio analysis in high dimensions: Residual-based nodewise regression in factor models. Journal of Econometrics 235 (2), 393-417.

Fan, J., Y. Liao, and M. Mincheva (2013). Large covariance estimation by thresholding principal orthogonal complements. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 75 (4), 603–680.

Ledoit, O, M. and M. Wolf (2017). Nonlinear shrinkage of the covariance matrix for portfolio selection: Markowitz meets goldilocks. Review of Financial Studies 30, 4349–4388.
