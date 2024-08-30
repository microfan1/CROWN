function w = GetWeights_WC(inv_Sigma, mu, wx, Delta, oneR)
    %inv_Sigma is the estimated precision matrix
    %mu is the mean of return
    %Delta is the aversion parameter of risk (see optimal weights under only weight constraints)
    %oneR is the selected stocks list
    %omega is the weights sum of selected stocksOne = ones(size(mu));
    One = ones(size(mu));
    a=inv_Sigma*One/(One'*inv_Sigma*One);
    B1 = One'*inv_Sigma*mu;
    kw = One'*inv_Sigma*mu/Delta;
    k = inv_Sigma*oneR/(One'*inv_Sigma*oneR);
    wk = oneR'*inv_Sigma*oneR/(oneR'*inv_Sigma*One);
    wa = oneR'*inv_Sigma*One/(One'*inv_Sigma*One);
    L = (k-a)/(wk - wa);
    wu = oneR'*inv_Sigma*mu/(One'*inv_Sigma*mu) - wa;
    w = (kw*(inv_Sigma*mu/B1-a)+(wx-kw*wu)*L) + (a-L*wa);
end
