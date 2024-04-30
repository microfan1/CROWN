# Nodewise estimation of the covariance matrix
est_ndwcov <- function(Y,ic){
  
  # initialization
  p <- ncol(Y)
  n <- nrow(Y)
  Y <- Y- t((apply(Y,2,mean))%*%matrix(1,1,n)) # Y is de-meaned
  
  C <- matrix(0,p,p)
  diag(C) <- 1
  tau <- NULL
  
  # Loop over the assets
  for(j in 1:p){
    # Estimate the Lasso
    jlas <- glmnet(x=Y[,-j],y=Y[,j],family = 'gaussian',intercept = FALSE)
    # Get fit
    jfit <- predict(jlas, newx=Y[,-j], type="response")    
    # residuals
    jres <- matrix(Y[,j],n,length(jlas$lambda)) - jfit
    # std err
    jsig <- colSums(jres^2)/n
    # Computing information criterion
    if(ic=='WIC') jbic  <- log(jsig) + jlas$df *log(n)/n * log(log(p)) # BIC (Wang,2010)
    if(ic=='GIC') jbic  <- log(jsig) + jlas$df *log(p)/n * log(log(n)) # GIC
    if(ic=='BIC') jbic  <- log(jsig) + jlas$df *log(n)/n  #BIC
    if(ic=='MIC') jbic  <- jsig + jlas$df *log(n) * log(log(p))/n # MC's IC
    if(ic=='AIC') jbic  <- log(jsig) + 2 * jlas$df # AIC
    # Index of selected model 
    jind  <- which.min(jbic)
    # Get the parameters
    jpar <- jlas$beta[,jind]
    # Computing tau squared
    jtau <- sum(jres[,jind]^2)/n + (1/2)*jlas$lambda[jind]*sum(abs(jpar)) # using (10)
    # Storing the parameters
    C[j,-j] <- -jpar
    tau <- c(tau,jtau)
  }
  
  # Construct T-squared inverse
  T2inv <- diag(1/tau)
  
  # Construct Theta-hat
  Theta <- T2inv %*% C
  
  # sparsity
  sp <- sum(Theta==0)/(p^2)
  
  return(list(NULL,Theta,sp))
}



est_ndwcov_factor <- function(Y,factors,ic,lambda.min = TRUE){
  
  # initialization
  p <- ncol(Y)
  n <- nrow(Y)
  C <- matrix(0,p,p)
  diag(C) <- 1
  tau <- NULL
  ns1 = rep(1,n)
  factormodel = lm(Y~factors)
  u = factormodel$residuals
  beta = t(coef(factormodel)[-1,])
  #Y = scale(Y,scale = FALSE)
  
  # Loop over the assets
  for(j in 1:p){
    
    if (ic!="cv"){
      # Estimate the Lasso
      jlas <- glmnet(x=u[,-j],y=u[,j],family = 'gaussian',standardize = FALSE,intercept = FALSE)
      # Get fit
      jfit <- predict(jlas, newx=u[,-j], type="response")    
      # residuals
      jres <- matrix(u[,j],n,length(jlas$lambda)) - jfit
      # std err
      jsig <- colSums(jres^2)/n
      # Computing information criterion
      if(ic=='WIC') jbic  <- log(jsig) + jlas$df * log(n)/n * log(log(p)) # BIC (Wang,2010)
      if(ic=='BIC') jbic  <- log(jsig) + jlas$df * log(n)/n  #BIC
      if(ic=='GIC') jbic  <- log(jsig) + jlas$df * log(p) * log(log(n))/n # Fan & Tang JRSS-B 2004
      if(ic=='AIC') jbic  <- log(jsig) + 2 * jlas$df # AIC 
      # Index of selected model 
      jind  <- which.min(jbic)
      jpar <- jlas$beta[,jind]
      jtau <- sum(u[,j]*jres[,jind])/n# + jlas$lambda[jind]*sum(abs(jpar)) # using (12)
    }else{
      jlas = cv.glmnet(x=u[,-j],y=u[,j],family = 'gaussian',standardize = FALSE,intercept = FALSE,nfolds = 5)
      if(lambda.min==TRUE){
        jfit <- predict(jlas, newx=u[,-j], lambda = "lambda.min")    
        jres <- u[,j] - jfit
        jpar = coef(jlas,lambda = "lambda.min")[-1]
        jtau = sum(u[,j]*jres)/n# + jlas$lambda.min*sum(abs(jpar))
      }else{
        jfit <- predict(jlas, newx=u[,-j], lambda = "lambda.1se")    
        jres <- u[,j] - jfit
        jpar = coef(jlas,lambda = "lambda.1se")[-1]
        jtau = sum(u[,j]*jres)/n# + jlas$lambda.min*sum(abs(jpar))
      }
      
    }
    
    C[j,-j] <- -jpar/jtau
    tau <- c(tau,jtau)
  }
  
  
  diag(C) = 1/tau
  omega = C
  omegasym = (C+t(C))/2
  
  covft = (1/n)*t(factors)%*%factors-(1/(n^2))*t(factors)%*%ns1%*%t(ns1)%*%factors 
  
  if(ncol(factors)==1){
    p1 = solve(solve(covft)+beta%*%omegasym%*%t(beta))
    TAU = omega - omega%*%t(beta)%*%p1%*%beta%*%omega
  }else{
    p1 = solve(solve(covft)+t(beta)%*%omegasym%*%beta)
    TAU = omega - omega%*%beta%*%p1%*%t(beta)%*%omega
  }
  
  
  return(TAU)
}

