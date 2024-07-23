function w = GetWeights_TE(inv_Sigma, mu, Xi, m)
    %inv_Sigma is the estimated precision matrix
    %mu is the mean of return
    %Xi is risk aversion factor
    One = ones(size(mu));
    a=inv_Sigma*One/(One'*inv_Sigma*One);
    t=inv_Sigma*mu/(One'*inv_Sigma*mu);
    u=t-a;
    kappa=One'*inv_Sigma*mu/Xi;
    w=kappa*u + m;
end
