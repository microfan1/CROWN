library(R.matlab)
library(glmnet)

#Remember to change this path
base_path <- 'E:/elective/RA/CUHK/Tutorial/code/methods/NODEWISE'

source(file.path(base_path, 'Nodewise_source.R'))

Data <- readMat(file.path(base_path, 'to_Nodewise_R.mat'))
y  =  Data$R

sigma_nw_list = est_ndwcov(y,'GIC')
sigma_nw = sigma_nw_list[[2]]

write.table(sigma_nw,
            file = file.path(base_path, 
                             'inv_Sigma_nw.csv'),
            col.names=F,
            row.names=F,
            sep=',') 

fac = Data$fac
sigma_nw_factor = est_ndwcov_factor(y,t(fac),'GIC')

write.table(sigma_nw_factor,
            file = file.path(base_path, 
                             'inv_Sigma_nw_factor.csv'),
            col.names=F,
            row.names=F,
            sep=',') 

