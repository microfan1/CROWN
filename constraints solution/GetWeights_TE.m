function w = GetWeights_TE(inv_Sigma, mu, Xi, m)
    %inv_Sigma is the estimated precision matrix
    %mu is the mean of return
    %Xi is risk aversion factor
    % m is the targeted index, e.g., equal/mv weight, known indices s&p 500, Russell 2000, etc.
    One = ones(size(mu));
    a=inv_Sigma*One/(One'*inv_Sigma*One);
    t=inv_Sigma*mu/(One'*inv_Sigma*mu);
    u=t-a;
    kappa=One'*inv_Sigma*mu/Xi;
    w=kappa*u + m;
end
