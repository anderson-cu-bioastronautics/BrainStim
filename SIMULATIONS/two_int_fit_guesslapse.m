function J = two_int_fit_guesslapse(x, X, cor, lambda_max, gamma_max)

lambda_min = 0;
% lambda_max = 0.1;
sigma_min = 0;

lambda_hat = x(1);
gamma_hat = x(2);
mu_hat = x(3);
sigma_hat = x(4);


% if gamma_hat < 0
%     gamma_hat = 0;
% end
% if lambda_hat < 0
%     lambda_hat = 0;
% end
% if sigma_hat < 0
%     sigma_hat = 0;
% end

F = cdf('norm', X, mu_hat, sigma_hat);
pred = 0.5 + gamma_hat + (1-0.5-gamma_hat-lambda_hat)*F;
% J = sum( (Y - pred).^2 );
J = -(sum(log(pred(cor==1))) + sum(log(1-pred(cor==0))));

% forceably apply limits to parameters by setting cost function to very
% large value
Jmax = length(X);
if lambda_hat<lambda_min || lambda_hat>lambda_max || sigma_hat<sigma_min || abs(gamma_hat)>gamma_max
%     display('Parameter out of range')
    J = Jmax;
end
    
    
    
    