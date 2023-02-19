%% T-test with estimated thresholds and Satterthwaite approximation
function[t_test] = t_test_satterthwaite(stim1,stim2,cor1,cor2,guesses,lambda_max)
% This function takes two sets of independently collected data used to
% estimate perception thresholds and estimates the p-value with a t-test.
% The t-test requires the use of the Sattherthwaite approximation for the
% degrees of freedom.

% Inputs:
%     stim1:      The stimulus values for the first set of data (in SR context, the sham set)
%     stim2:      The stimulus values for the second set of data (in SR context, the optimal set)
%     cor1:       The vector of responses for the first set of data
%     cor2:       The vector of responses for the second set of data
%     guesses:    A 1x3 vector containing initial guesses for fitting the psychometric curves
%     lambda_max: The maximum allowable value of lambda, used for fitting the psychometric curves
%     
% Outputs:
%     x1:         The fitted threshold for set 1
%     x2:         The fitted threshold for set 2
%     se1:        The standard error for set 1
%     se2:        The standard error for set 2
%     dof1:       Degrees of freedom for set 1
%     dof2:       Degrees of freedom for set 2
%     dof:        Overall degrees of freedom
%     t:          The calculated t-value
%     p:          The calculated p-value

%% Estimate the thresholds

% "Sham" or base set (set 1)
params1 = fminsearch(@(params1) two_int_fit_simp(params1, stim1, cor1), guesses(1:2));

% Extract the fitted threshold
x1 = params1(1);

% Comparison set (set 2)
params2 = fminsearch(@(params2) two_int_fit_simp(params2, stim2, cor2), guesses(1:2));

% Extract the fitted threshold
x2 = params2(1);

%% Calculate the standard errors using Jackknife error estimation

% How many trials did we have?
N1 = length(stim1);
N2 = length(stim2);

% Initialize vectors to contain the different results of mu for set 1
% jackknife
mu1 = zeros(1,N1);

% Use the jackknife method to estimate the threshold for each trial removed
for i = 1:N1
    
    % Make a copy of the original set
    stim1_copy = stim1;
    cor1_copy = cor1;
    
    % Remove the ith trial
    stim1_copy(i) = [];
    cor1_copy(i) = [];
    
    % Generate the fit with the ith trial removed
    paramsi = fminsearch(@(paramsi) two_int_fit_simp(paramsi, stim1_copy, cor1_copy), guesses(1:2));
    
    % Fill in the vector
    mu1(i) = paramsi(2);
end

% Calculate the standard error for set 1
% se1 = sqrt((N1-1)/N1*sum( (mu1 - x1*ones(size(mu1))).^2));
se1 = sqrt(sum( (mu1 - x1*ones(size(mu1))).^2)/(N1*(N1-1)));

% Initialize vector to contain the different results fo mu for set 2
% jackknife
mu2 = zeros(1,N2);

% Use the jackknife method to estimate the threshold for each trial removed
for i = 1:N2
    
    % Make a copy of the original set
    stim2_copy = stim2;
    cor2_copy = cor2;
    
    % Remove the ith trial
    stim2_copy(i) = [];
    cor2_copy(i) = [];
    
    % Generate the fit with the ith trial removed
    paramsi = fminsearch(@(paramsi) two_int_fit_simp(paramsi, stim2_copy, cor2_copy), guesses(1:2));
    
    % Fill in the vector
    mu2(i) = paramsi(2);
end

% Calculate the standard error for set 2
temp1 = mu2 - x2*ones(size(mu2));
temp2 = temp1.^2;
temp3 = sum(temp2);
se2 = sqrt(sum( (mu2 - x2*ones(size(mu2))).^2)/(N2*(N2-1)));

%% Calculate the t-value to use in the t-test

% Calculate the t-value
t = abs((x1 - x2)/sqrt((se1^2)/N1 + (se2^2)/N2));

%% Calculate the DOF for each set of data using the Satterthwaite approximation

% Difference between overall threshold and each new estimate from jackknife
mu_diff1 = ones(1,length(mu1))*x1 - mu1;
mu_diff2 = ones(1,length(mu2))*x2 - mu2;

% The numerator terms
numerator1 = (sum(mu_diff1.^2)^2);
numerator2 = (sum(mu_diff2.^2)^2);

% The denominator terms
denominator1 = sum(mu_diff1.^4);
denominator2 = sum(mu_diff2.^4);

% Each DOF
dof1 = numerator1/denominator1;
dof2 = numerator2/denominator2;

%% Calculate the overall DOF

% The numerator for the overall DOF
numerator = ((se1^2 + se2^2)^2);

% The demonimator for the overall DOF
denominator = (se1^4/dof1) + (se2^4/dof2);

% The overall DOF
dof = numerator/denominator;

%% Get the p-value

p = 1-tcdf(t,dof);

t_test = struct();
t_test.X1 = x1;
t_test.X2 = x2;
t_test.SE1 = se1;
t_test.SE2 = se2;
t_test.DOF1 = dof1;
t_test.DOF2 = dof2;
t_test.DOF = dof;
t_test.t = t;
t_test.p = p;

end