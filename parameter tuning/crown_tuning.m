% subdata  is N by p matrix
% this function return Nodewise regression Precision matrix estimator (Crown)
function [Pm,C_hat,T_2_hat] = crown_tuning(subdata, tuning_method)

if nargin > 1
   valid_methods = {'GIC', 'BIC', 'AIC', 'EBIC', 'CV'};
   if ~ismember(tuning_method, valid_methods)
        error('Unsupported tuning method. Use one of ''GIC'', ''BIC'', ''AIC'', ''EBIC'' or ''CV''.');
   end
else
   tuning_method = 'GIC';
end

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
    
    if strcmp(tuning_method, 'CV')
        K = 5; 
        for k = 1:K
            indices = randperm(n);
            val_begin = (k-1)*floor(n/K)+1;
            val_end = val_begin + floor(n/K)-1;
            val_indices = [];
            val_indices = indices(val_begin:val_end);
            train_indices = [];
            train_indices = setdiff(indices, val_indices);
            fit_temp = glmnet(X(train_indices,:), Y(train_indices), [], option1);
            gamma_temp = (fit_temp.beta)';
            lambda_set = fit_temp.lambda;
            mse_cv = zeros(size(lambda_set));
            for i = 1:length(lambda_set)
                gamma = (gamma_temp(i,:))';
                e_temp = Y(val_indices) - X(val_indices,:) * gamma;
                mse_cv(i) = e_temp' * e_temp / length(val_indices); 
            end
        end
        mse_cv = mse_cv / K;
        [~, aim_loc] = min(mse_cv);
        aim_fit = gamma_temp(aim_loc,:);
        aim_lambda = lambda_set(aim_loc);
    else
        fit2 = glmnet(X,Y,[],option1);  
        fit1 = (fit2.beta)';
        lambda_set = fit2.lambda;
        IC_set = [];
        sigma_set = [];
        for ic_loop = 1 : length(lambda_set)
            gamma = (fit1(ic_loop,:))';
            e = Y_0 - X_0* gamma;
            sigma_2 = e'*e/n;
            S_abs = sum(gamma ~= 0);
            sigma_set(ic_loop) = sigma_2;
            switch tuning_method
                case 'GIC'
                    GIC = log(sigma_2) + S_abs * (log(p) * log(log(n))) / n;
                    IC_set(ic_loop) = GIC;
                case 'BIC'
                    BIC = log(sigma_2) - (S_abs * log(n)) / 2;
                    IC_set(ic_loop) = BIC;
                case 'AIC'
                    AIC = log(sigma_2) + (2 * S_abs) / n;
                    IC_set(ic_loop) = AIC;
                case 'EBIC'
                    EBIC = log(sigma_2) + S_abs * (log(n) + 0.2 * log(p)) / n;
                    IC_set(ic_loop) = EBIC;
            end
        end
        [~,aim_loc] = min(IC_set);
        aim_fit = fit1(aim_loc,:);
        aim_lambda = lambda_set(aim_loc);
    end

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