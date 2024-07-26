function [w_err,te,risk_hat,V_err,SR_hat,SRerr]=GetPerformance(w_test,w_pop,V,R,mu,m)  
    % w_test is the weights that you wanna test performance
    % w_pop is the populational solution
    % V is portfolio variance
    % R is return matrix
    % mu is the mean of the return
    One = ones(size(mu));
    %Weight Error
    w_err = norm(w_test-w_pop,1);
    % Test empirical version, we should use weight itself
    te =  std(R*w_test - R*m);
    %\hat risk
    risk_hat   =  sqrt(w_test'*V*w_test);
    % Covariance error
    risk_real =  sqrt(w_pop'*V*w_pop);
    V_err    =  abs((risk_hat/risk_real)^2 -1);
    %\hat Sharpe ratio
    SR_real  =  w_pop'*mu/risk_real;
    SR_hat = w_test'*mu/risk_hat;
    %Sharpe Ratio Error
    SRerr =  abs((SR_hat/SR_real)^2 - 1);
end
