function w = GetWeights_NC(inv_Sigma, mu, delta)
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
