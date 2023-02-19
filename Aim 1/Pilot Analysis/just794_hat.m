% function handle for fitting lapse model

function [musi_hat]=just794_hat(X,cor,mu_guess,sigma_guess)

guess_rate = 0.5;  % the level at which subject guesses correctly, normally 0.5

lambda_max = 0.06;
% gamma_max = 0.06;
x = fminsearch(@(x) two_int_fit_lapse(x, X, cor, lambda_max), [0 mu_guess sigma_guess]);
lambda_hat = x(1);
mu_hat = x(2);
sigma_hat = x(3);
musi_hat=[mu_hat sigma_hat lambda_hat];
% thresh_794_hat_boot = icdf('norm', (0.794 - guess_rate)/(1-0.5-lambda_hat), mu_hat, sigma_hat);

