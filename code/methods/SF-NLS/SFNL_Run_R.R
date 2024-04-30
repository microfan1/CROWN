library(nlshrink)
library(R.matlab)
library(PDSCE)

#Remember to change this path
base_path <- 'E:/elective/RA/CUHK/Tutorial/code/methods/SF-NLS'
Data <- readMat(file.path(base_path, 'to_SF-NLS_R.mat'))
Y <- Data$R

fact = rowMeans(Y)
#fact = t(Data$fac)

factm = lm(Y ~ fact)
loadings = matrix(coef(factm)[2,])
res = residuals(factm)
variances = apply(res,2,var)

variances = diag(variances)

sigmaf = loadings%*%var(fact)%*%t(loadings) + variances
sqrtsigf = expm::sqrtm(sigmaf)
sqrtsigfinv = chol2inv(chol(sqrtsigf))

Yadj = Y%*% sqrtsigfinv

sigmahat = nlshrink_cov(Yadj, k=1) # compute non-linear shrinkage estimate

sigma = sqrtsigf %*% sigmahat%*%sqrtsigf


write.table(sigma,
            file = file.path(base_path, 
                             'Sigma_SF-NLS.csv'),
            col.names=F,row.names=F,sep=',') 
