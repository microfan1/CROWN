HT4SR_bs <- function(returns, benchmark, batch=100, nBoot=500, bBoot=1){
  #' @details A function performs the testing of Sharpe ratio difference for 
  #' multiple methods (in our case, NW, CROWN, POET, NLS, and SFNL) using bootstrap
  #' (https://cran.r-project.org/web/packages/PeerPerformance/index.html).
  #' Since the `sharpeScreening` function is performing a t-test on the difference
  #' of two portfolios, dividing it by 2 yields a p-value for a one-side test
  #' where
  #'  H0: SR_1 - SR_2 <= 0; H1: SR_1 - SR_2 > 0.
  #'  That is, we are interesting in whether the SR of our benchmark
  #'  portfolio is significantly higher than that of the counterpart.
  #'  
  #' @param returns A T-by-N matrix of T returns for N funds
  #' @param benchmark The column of the benchmark portfolio (in our case, CROWN)
  #' @param batch Number of iterations of bootstrap with nBoot
  #' @param nBoot Number of bootstrap replications, default=500
  #' @param bBoot Block length, default = 1
  #' @return  A vector of average p-values of test of Sharpe ratio differences
  #' across batches
  
  p_list = matrix(NA, nrow = batch, ncol = ncol(returns))
  for(i in 1:batch){
    ctr = list(type = 2, nBoot = nBoot, bBoot = bBoot)
    p_matrix = sharpeScreening(as.matrix(returns), control = ctr)$pval
    p_list[i, ] = p_matrix[benchmark, ]
  }
  output = colMeans(p_list)
  return(output / 2)
}