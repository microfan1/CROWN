function w = GetWeights_TEWC(inv_Sigma, mu, Xi, oneR, omega)
    %inv_Sigma is the estimated precision matrix
    %mu is the mean of return
    %Xi is risk aversion factor
    %oneR is the selected stocks list
    %omega is the weights sum of selected stocks
    One = ones(size(mu));
    a=inv_Sigma*One/(One'*inv_Sigma*One);
    t=inv_Sigma*mu/(One'*inv_Sigma*mu);
    u=t-a;
    kappa=One'*inv_Sigma*mu/Xi;
    w0=kappa*u;
    k = inv_Sigma*oneR/(One'*inv_Sigma*oneR);
    wk = oneR'*inv_Sigma*oneR/(oneR'*inv_Sigma*One);
    wa = oneR'*inv_Sigma*One/(One'*inv_Sigma*One);
    L = (k-a)/(wk - wa);
    wu = oneR'*inv_Sigma*mu/(One'*inv_Sigma*mu) - wa;
    w = (omega - kappa*wu)*L + w0;
end
