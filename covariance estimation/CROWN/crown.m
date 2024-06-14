
function [Pm,C_hat,T_2_hat] = crown(subdata)
% subdata  is N by p matrix
% this function return Nodewise regression Precision matrix estimator(Crown)
n = size(subdata,1);
p = size(subdata,2);
r_bar = mean(subdata);
subdata_star = subdata - repmat(r_bar,n,1);
aim_fit_set = [];
T_2 = [];
for b_loop = 1 : p
    Y = subdata_star(:,b_loop);
    X = subdata_star;
    X(:,b_loop) = [];
    Y_0 = subdata(:,b_loop);
    X_0 = subdata;
    X_0(:,b_loop) = [];
    option1 = glmnetSet;
    option1.intr = 0;
    option1.standardize = 0;
    fit2 = glmnet(X,Y,[],option1);  
    fit1 = (fit2.beta)';
    lambda_set = fit2.lambda;
    GIC_set = [];
    sigma_set = [];
    for gic_loop = 1 : length(lambda_set)
        gamma = (fit1(gic_loop,:))';
        e = Y_0 - X_0* gamma;
        sigma_2 = e'*e/n;
        S_abs = sum(gamma ~= 0);
        GIC = log(sigma_2) + S_abs * (log(p)* log(log(n)))/ n;
        sigma_set(gic_loop) = sigma_2;
        GIC_set(gic_loop) = GIC;
    end
    [~,aim_loc] = min(GIC_set);
    aim_fit = fit1(aim_loc,:);
    aim_lambda = lambda_set(aim_loc);
    e_star = Y - X * aim_fit';
    tau = (Y' * e_star)/n; 
    T_2(b_loop) = tau;
    aim_fit_set(b_loop,:) = aim_fit;
end
C_hat = eye(p);
for i = 1:p
    for j = 1:p
        if i == j
            continue;
        end
        if i < j
            C_hat(i,j) = -aim_fit_set(i,j-1);
        end
        if i > j
            C_hat(i,j) = -aim_fit_set(i,j);
        end
    end
end
T_2_hat = diag(T_2);
Pm = (T_2_hat) \ C_hat;
for i = 1 : size(Pm,1)
    for j = i : size(Pm,1)
        if abs(Pm(i,j)) < abs(Pm(j,i))
           Pm(j,i) = Pm(i,j);
        else
            Pm(i,j) = Pm(j,i);
        end
    end
end
end