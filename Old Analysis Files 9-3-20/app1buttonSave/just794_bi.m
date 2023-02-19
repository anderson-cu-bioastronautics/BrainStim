% function handle for fitting binary model

function [musi_bi]=just794_bi(X,cor,mu_guess,sigma_guess)
guess_rate = 0.5;  % the level at which subject guesses correctly, normally 0.5

x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mu_guess, sigma_guess]);
mu_bi = x(1);
sigma_bi = x(2);
musi_bi=[mu_bi sigma_bi];
% thresh_794_bi_boot = icdf('norm', (0.794 - guess_rate)/(1-guess_rate), mu_bi, sigma_bi);
