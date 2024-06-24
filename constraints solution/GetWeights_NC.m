function w = GetWeights_NC(inv_Sigma, mu, delta)
    %It's worth mentioning that the weights we get here is an adjusted version, which means it should be added by tracking index b then becomes the real one.
    %This function is used to get weights under no constraints 
    %inv_Sigma is the estimated precision matrix
    %mu is the mean of return
    %delta is the aversion parameter of risk
    One = ones(size(mu));
    a=inv_Sigma*One/(One'*inv_Sigma*One);
    B1 = One'*inv_Sigma*mu;
    kw = One'*inv_Sigma*mu/delta;
    w = kw*(inv_Sigma*mu/B1-a) + a;
end
