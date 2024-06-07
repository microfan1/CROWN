function w = GetWeights_TEOnly(inv_Sigma, mu, Xi)
    %inv_Sigma is the estimated precision matrix
    %mu is the mean of return
    %Xi is risk aversion factor
    One = ones(size(mu));
    a=inv_Sigma*One/(One'*inv_Sigma*One);
    t=inv_Sigma*mu/(One'*inv_Sigma*mu);
    u=t-a;
    kappa=One'*inv_Sigma*mu/Xi;
    w=kappa*u;
end