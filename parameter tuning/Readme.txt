% Users can get the CROWN method (residual nodewise) precision matrix estimate with this function crown_tuning, which gives users the option to choose the tuning method including cross validations and several information criterion. The default tuning method is 'GIC'.

%Example: set the tuning method to (5 fold) cross validation
O=crown_tuning(res, 'CV');
inv_Sigma_crown=O-O*beta/(inv(cov_f)+beta'*(O+O')/2*beta)*beta'*O;
