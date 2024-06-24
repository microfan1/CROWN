function w = GetWeights_TE(inv_Sigma, mu, Xi)
    %It's worth mentioning that the weights we get here is an adjusted version, which means it should be added by tracking index b then becomes the real one.
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
