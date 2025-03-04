CROWN
================

## 1. Overview
This repository provides the code needed to compute the **C**onstrained **R**esidual Nodewise **O**ptimal **W**eight Regressio**n** (CROWN), which can be used to construct constrained portfolios in a high-dimensional scenario (p>T). It also shows different state-of-the-art methods in obtaining a solution for high dimensional portfolios with a large number of assets. Regarding the constraints, specifically, we consider portfolios with tracking error constraints, portfolios with tracking error jointly with weight (equality or inequality) restrictions, and portfolios with only weight restrictions. 

Firstly, we need an estimator of covariance/precision matrix for the return data.

This repository gives the implementation code for **CROWN** which uses a residual-based nodewise regression (Caner, Medeiros and Vasconcelos, 2023, J. Econom.) to get the estimate of the covariance. For comparison it also includes 4 other popular methods:
1. **nodewise** from Callot et al., 2021, [A Nodewise Regression Approach to Estimating Large Portfolios]
2. **POET** from Fan et al., 2013, [Large Covariance Estimation by Thresholding Principal Orthogonal Complements]
3. **NLS**  and 4. **SF-NLS** from Ledoit and Wolf, 2017, [Nonlinear Shrinkage of the Covariance Matrix for Portfolio Selection: Markowitz Meets Goldilocks]

The paper detailing the methodology for the CROWN estimator is available at: [Portfolio Analysis in High Dimensions with Tracking Error and Weight Constraints](https://arxiv.org/abs/2402.17523).

The S&P 500 stocks return data is from CRSP and accessable upon subscription to the WRDS. Data on historical index components is from Refinitiv.

Data dictionary:
https://wrds-www.wharton.upenn.edu/pages/get-data/center-research-security-prices-crsp/annual-update/stock-security-files/monthly-stock-file/

https://wrds-www.wharton.upenn.edu/pages/get-data/center-research-security-prices-crsp/annual-update/stock-security-files/stock-market-indexes/

https://wrds-www.wharton.upenn.edu/pages/get-data/fama-french-portfolios-and-factors/fama-french-portfolios/factors-monthly-frequency/


## 2. Required packages
You need to download matlab toolbox RunRCode for R functions calling.

Wei-Rong Chen (2024). RunRcode(RscriptFileName,Rpath) (https://www.mathworks.com/matlabcentral/fileexchange/50071-runrcode-rscriptfilename-rpath), MATLAB Central File Exchange. Retrieved July 24, 2024.

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
This repository contains the Matlab and R scripts used for the implementation of the aforementioned 5 methods. 

- covariance estimation/
   - CROWN: 
     - glmnet_matlab: Package used by `crown.m`.
     - `crown.m`: Function to get crown estimator.
   - nodewise: 
     - `Nodewise_Run_R.R`: Function to get Nodewise estimator.
     - `Nodewise_source.R`: Functions used by `Nodewise_Run_R.R`.
   - POET:
     - POET: Package used by `POET_R.R`.
     - `POET.m`: Code to link R and Matlab.
     - `POET_R.R`: Code to get POET estimator.
   - NLS:
     - `nls_covMarket.m`: Function to get nonlinear shrinkage estimator.
   - SF-NLS:
     - `SF-NLS_Run_R.R`: Code to get single-factor nonlinear shrinkage estimator.

## 5. Examples of Portfolio Construction under Constraints
`Example.mlx` exemplifies how to use CROWN and four other baseline methods for constructing portfolios with different kinds of constraints. Specifically, we consider portfolios with tracking error (TE) constraints, portfolios with tracking error jointly with weight constraints (TEWC), and portfolios with only weight constraints (WC). 

Please modify the working path and Rpath in `Example.mlx`, `Nodewise_Run_R.R`, `POET_R.R` and `SF-NLS_Run_R.R` first.

### 5.1 Generation of Example Data: 
1. 180 stocks (assets) with 150 trading records;
2. 3 factors;
3. Set Tracking Error constraint (e.g., TE<=0.5);
4. Set weight constraint, e.g., sum of first 10 stocks weights equals to 0.2 (i.e., the restricted vector of assets takes value of 1 for the first 10 assets, and 0 for the rest 170 assets, and wx=0.2).
5. Set delta equals to 0.5 for risk aversion parameter in weight-constraint-only problem.


For calibrated data simulations, just change the parameters (or load the real data in similar format):

```
% Set the portfolio and sample size
pNpairs = [180,150];p=pNpairs(1);N=pNpairs(2);

% Set the Tracking Error Constraint (e.g., required by the fund manager)
TEC = 0.5;

% Set the Weights Constraint, here we restrict the first 10 stocks, whose weights sum up to 0.2
oneR = [ones(10,1);zeros(p-10,1)];
wx = 0.2;

% Define Xi the risk aversion parameter for the given tracking error, in practice, Xi is a user defined value.
k0 = (One'*inv_Sigma_pop*mu_pop)/(One'*inv_Sigma_pop*One);
k0error = mu_pop-k0*One;
Xi = sqrt(k0error'*inv_Sigma_pop*k0error)/TEC;

% Define delta as risk aversion when only consider weights constraint, in practice, delta is a user defined value.
delta = 0.5;

% 3-Factors
flag=3; k=3;
% The following data loading can be replaced by user's own data, such as FF-factors or Q factors.
% factor loadings
mu_b=[-0.1,0.1,0.1];%kx1
factor_loading_b=mvnrnd(mu_b,0.1.*eye(3),p);
% the factors follow AR(1) process
alpha1 = 0.03;alpha2 = -0.05;alpha3 = -0.05;
alphamtx = [alpha1,0,0;0,alpha2,0;0,0,alpha3];% alpha matrix
mu_f = [0.005,0.005,0.005]';
cov_f = eye(3);
```

Use `Data_Generation.m` to generate simulation data based on these parameters. 

### 5.2 Estimation of Precision Matrix
The code to estimate precision matrix by different methods:
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

### 5.3 Get Weights and Report Performance

1. Only Consider Tracking Error Constraint
   
| Method | TE       | Weight_ER | Risk_ER    | SR       | SR_ER    |
|--------|----------|-----------|------------|----------|----------|
| NW     | 6.5155   | 16.0782   | 136.2974   | **0.0084**   | 0.1522   |
| CROWN  | **0.4932**   | **0.3853**    | **0.00019475** | 0.0078   | **0.0045**   |
| POET   | 1.9888   | 13.8279   | 22.1934    | 0.0074   | 0.1164   |
| NLS    | 0.8292   | 9.4759    | 9.4595     | 0.0074   | 0.1198   |
| SF-NLS | 0.8153   | 6.2999    | 5.2599     | 0.0077   | 0.0502   |

2. Only consider Weights Constraint


| Method | TE       | Weight_ER | Risk_ER   | SR      | SR_ER    |
|--------|----------|-----------|-----------|---------|----------|
| NW     | 0.2081   | 0.7663    | 5.4120    | **0.0064**  | 25.9909  |
| CROWN  | **0.2020**   | **0.2402**    | **0.0453**    | 0.0012  | **0.0228**   |
| POET   | 0.2209   | 0.7259    | 0.5171    | 0.0032  | 5.6544   |
| NLS    | 0.2160   | 0.8769    | 0.6057    | 0.0026  | 3.3970   |
| SF-NLS | 0.2263   | 0.4601    | 0.2181    | 0.0022  | 2.0982   |

3. Tracking Error + Weights Constraint

| Method | TE       | Weight_ER | Risk_ER   | SR      | SR_ER    |
|--------|----------|-----------|-----------|---------|----------|
| NW     | 6.4798   | 16.0437   | 133.1197  | **0.0084**  | 0.1678   |
| CROWN  | **0.4964**   | **0.4223**    | **0.0012**    | 0.0078  | **0.0052**   |
| POET   | 1.9786   | 13.8588   | 21.8348   | 0.0074  | 0.1048   |
| NLS    | 0.8287   | 9.4669    | 9.3650    | 0.0074  | 0.1097   |
| SF-NLS | 0.8166   | 6.3013    | 5.2000    | 0.0076  | 0.0400   |


4. Tracking Error +  Weights Constraint though non-binding 

| Method | TE       | Weight_ER | Risk_ER    | SR       | SR_ER    |
|--------|----------|-----------|------------|----------|----------|
| NW     | 6.5155   | 16.0782   | 136.2974   | **0.0084**   | 0.1522   |
| CROWN  | **0.4932**   | **0.3853**    | **0.00019475** | 0.0078   | **0.0045**   |
| POET   | 1.9888   | 13.8279   | 22.1934    | 0.0074   | 0.1164   |
| NLS    | 0.8292   | 9.4759    | 9.4595     | 0.0074   | 0.1198   |
| SF-NLS | 0.8153   | 6.2999    | 5.2599     | 0.0077   | 0.0502   |

## References

Callot, L., M. Caner, O. Onder, and E. Ulasan (2021). A nodewise regression approach to estimating large portfolios. Journal of Business and Economic Statistics 39, 520–531.

Caner, M., Fan, Q., and Li, Y. (2024). Navigating Complexity: Constrained Portfolio Analysis in High Dimensions with Tracking Error and Weight Constraints. arXiv preprint arXiv:2402.17523.

Caner, M., Medeiros, M., and G. Vasconcelos (2023). Sharpe Ratio analysis in high dimensions: Residual-based nodewise regression in factor models. Journal of Econometrics 235 (2), 393-417.

Fan, J., Y. Liao, and M. Mincheva (2013). Large covariance estimation by thresholding principal orthogonal complements. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 75 (4), 603–680.

Ledoit, O, M. and M. Wolf (2017). Nonlinear shrinkage of the covariance matrix for portfolio selection: Markowitz meets goldilocks. Review of Financial Studies 30, 4349–4388.
