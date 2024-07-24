You can replace the crown with this new function crown_tuning, which is equipped with parameter tuning by cross validation and several information criterion. The default tuning method is GIC.

Example:
O=crown_tuning(res, 'CV');
inv_Sigma_crown=O-O*beta/(inv(cov_f)+beta'*(O+O')/2*beta)*beta'*O;
