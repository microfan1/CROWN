function w = GetWeights_WC(inv_Sigma, mu, wx, delta, oneR)
    %inv_Sigma is the estimated precision matrix
    %mu is the mean of return
    %delta is the risk aversion parameter
    %oneR is the selected constrained stocks list
    %omega is the weights sum of selected stocks;
    One = ones(size(mu));
    a=inv_Sigma*One/(One'*inv_Sigma*One);
    B1 = One'*inv_Sigma*mu;
    kw = One'*inv_Sigma*mu/delta;
    k = inv_Sigma*oneR/(One'*inv_Sigma*oneR);
    wk = oneR'*inv_Sigma*oneR/(oneR'*inv_Sigma*One);
    wa = oneR'*inv_Sigma*One/(One'*inv_Sigma*One);
    L = (k-a)/(wk - wa);
    wu = oneR'*inv_Sigma*mu/(One'*inv_Sigma*mu) - wa;
    w = (kw*(inv_Sigma*mu/B1-a)+(wx-kw*wu)*L) + (a-L*wa);
end
