Crown
================

## 1. Overview
This repository provides the code needed to compute the CROWN method for **C**onstrained **R**esidual Nodewise **O**ptimal **W**eight Regressio**n**, which can be used to construct constrained portfolios in a high-dimensional scenario (p>T). It also shows different state-of-the-art methods in obtaining a solution for high dimensional portfolios with a large number of assets. Regarding the constraints, specifically, we consider portfolios with tracking error constraints, portfolios with tracking error jointly with weight (equality or inequality) restrictions, and portfolios with only weight restrictions. 

Firstly, we need an estimator of covariance/precision matrix for the return data.

This repository shows the method **CROWN** using a residual-based nodewise regression (Caner, Medeiros and Vasconcelos, 2023, J. Econom.) to get the estimate of the covariance, along with 4 other popular methods:
1. **nodewise** from [A Nodewise Regression Approach to Estimating Large Portfolios]
2. **POET** from [Large Covariance Estimation by Thresholding Principal Orthogonal Complements]
3. **NLS**  from [Nonlinear Shrinkage of the Covariance Matrix for Portfolio Selection: Markowitz Meets Goldilocks]
4. **SF-NLS** from [Nonlinear Shrinkage of the Covariance Matrix for Portfolio Selection: Markowitz Meets Goldilocks]

The paper detailing the methodology for the CROWN estimator is available at: [Navigating Complexity: Constrained Portfolio Analysis in High Dimensions with Tracking Error and Weight Constraints](https://arxiv.org/abs/2402.17523).


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
     - `nls_covMarket.m`: Function to get nonlinear shrinkage estimator.
   - nodewise: 
     - `Nodewise_Run_R.R`: Function to get Nodewise estimator.
     - `Nodewise_source.R`: Functions used by `Nodewise_Run_R.R`.
   - POET:
     - POET: Package used by `POET_R.R`.
     - `POET.m`: Code to link R and Matlab.
     - `POET_R.R`: Code to get POET estimator.
   - SF-NLS:
     - `SF-NLS_Run_R.R`: Code to get single-factor nonlinear shrinkage estimator.
## 5. Example of Usage
`Example.mlx` shows how to use CROWN and four other baseline methods for constructing portfolios with different kinds of constraints. Specifically, we consider portfolios with tracking error (TE) constraints, portfolios with tracking error jointly with weight constraints (TEWC), and portfolios with only weight constraints (WC). 

Please modify the working path and Rpath in `Example.mlx`, `Nodewise_Run_R.R`, `POET_R.R` and `SF-NLS_Run_R.R` first.

### 5.1 Example and Data Loading: 
1. 180 stocks with 150 trading records;
2. 3 factors;
3. Set Tracking Error constraint equals 0.5 (i.e., TE<=0.5);
4. Set sum of first 10 weights equals to 0.2 for weights constraint (i.e., the restricted vector of assets takes value of 1 for the first 10 assets, and 0 for the rest 170 assets, if the benchmark index has first 10 .
5. Set delta equals to 0.5 for risk aversion parameter in weight-constraint-only problem.


For data loading below, feel free to change them (or load your own data in similar format):

```
% Set your portfolio and sample size
pNpairs = [180,150];p=pNpairs(1);N=pNpairs(2);

% Set your Tracking Error Constraint
TEC = 0.5;

% Set your Weights Constraint
oneR = [ones(10,1);zeros(p-10,1)];
omega = 0.2;

% Define Xi for tracking error
k0 = (One'*inv_Sigma_pop*mu_pop)/(One'*inv_Sigma_pop*One);
k0error = mu_pop-k0*One;
Xi = sqrt(k0error'*inv_Sigma_pop*k0error)/TEC;

%Define delta as risk aversion when only consider weights constraint
delta = 0.5;

%3-Factor
flag=3;k = 3;
%This Data Loading is used in Crown Draft 
beta_total=mvnrnd([0.005,0.005,0.005],0.1.*eye(3),p);
bf = [0.03,0,0;0,-0.05,0;0,0,-0.05];% B Ft + et
mu_f=[-0.1,0.1,0.1]';%kx1
cov_f = eye(3); % Used for simulation
```

Use `Data_Generation.m` to generate simulation data based on these parameters. 

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
  O=crown(res);
  inv_Sigma_crown=O-O*beta/(inv(cov_f)+beta'*(O+O')/2*beta)*beta'*O;

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

### 5.3 Get Weights and Performance
It's worth mentioning that the TE constrained weights we get fron 'Get_weights' functions are the differences between the studied portfolio and a benchmark index, which means they should be added by tracking index b to get the real weights.
1. Only Consider Tracking Error Constraint
   
| Method | TE     | Weight_ER | Risk_ER   | SR      | SR_ER   |
|--------|--------|-----------|-----------|---------|---------|
| NW     | 5.9653 | 14.9595   | 125.8946  | **0.1678**  | 0.0185  |
| CROWN  | **0.5126** | **0.4497**    | **0.00004210**| 0.1658  | **0.0065**  |
| POET   | 1.2842 | 9.2783    | 13.7568   | 0.1432  | 0.2589  |
| NLS    | 0.7456 | 10.1683   | 11.3073   | 0.1435  | 0.2555  |
| SF-NLS | 0.7501 | 5.6867    | 5.2606    | 0.1510  | 0.1755  |

2. Only consider Weights Constraint


| Method | TE      | Weight_ER | Risk_ER  | SR      | SR_ER   |
|--------|---------|-----------|-----------|---------|---------|
| NW     | 4.0267  | 10.1440   | 129.7598  | **0.1678**  | 0.0867  |
| CROWN  | **0.3299**  | **0.3252**| **0.0014**| 0.1602  | **0.0092**  |
| POET   | 0.8269  | 6.0460    | 13.3640   | 0.1430  | 0.2105  |
| NLS    | 0.4802  | 6.9009    | 11.1970   | 0.1419  | 0.2232  |
| SF-NLS | 0.4741  | 3.8202    | 5.0338    | 0.1492  | 0.1404  |

3. Tracking Error + Weights Constraint

| Method | TE     | Weight_ER | Risk_ER   | SR      | SR_ER   |
|--------|--------|-----------|-----------|---------|---------|
| NW     | 5.9660 | 14.9389   | 123.1144  | **0.1678**  | 0.0396  |
| CROWN  | **0.5164** | **0.4726**    | **0.00005900**| 0.1640  | **0.0077**  |
| POET   | 1.2796 | 9.1145    | 13.2629   | 0.1434  | 0.2412  |
| NLS    | 0.7457 | 10.1618   | 11.0449   | 0.1435  | 0.2403  |
| SF-NLS | 0.7497 | 5.6423    | 5.1530    | 0.1509  | 0.1599  |


4. Tracking Error +  Weights Constraint though non-binding 

| Method | TE     | Weight_ER | Risk_ER   | SR      | SR_ER   |
|--------|--------|-----------|-----------|---------|---------|
| NW     | 5.9653 | 14.9595   | 125.8946  | **0.1678**  | 0.0185  |
| CROWN  | **0.5126** | **0.4497**    | **0.00004210**| 0.1658  | **0.0065**  |
| POET   | 1.2842 | 9.2783    | 13.7568   | 0.1432  | 0.2589  |
| NLS    | 0.7456 | 10.1683   | 11.3073   | 0.1435  | 0.2555  |
| SF-NLS | 0.7501 | 5.6867    | 5.2606    | 0.1510  | 0.1755  |

## References

Callot, L., M. Caner, O. Onder, and E. Ulasan (2021). A nodewise regression approach to estimating large portfolios. Journal of Business and Economic Statistics 39, 520–531.

Caner, M., Fan, Q., and Li, Y. (2024). Navigating Complexity: Constrained Portfolio Analysis in High Dimensions with Tracking Error and Weight Constraints. arXiv preprint arXiv:2402.17523.

Caner, M., Medeiros, M., and G. Vasconcelos (2023). Sharpe Ratio analysis in high dimensions: Residual-based nodewise regression in factor models. Journal of Econometrics 235 (2), 393-417.

Fan, J., Y. Liao, and M. Mincheva (2013). Large covariance estimation by thresholding principal orthogonal complements. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 75 (4), 603–680.

Ledoit, O, M. and M. Wolf (2017). Nonlinear shrinkage of the covariance matrix for portfolio selection: Markowitz meets goldilocks. Review of Financial Studies 30, 4349–4388.
