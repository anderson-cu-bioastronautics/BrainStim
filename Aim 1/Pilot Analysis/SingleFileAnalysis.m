% Program to read in and analyze one data file

% Fid is the filename to extract data from, must provide estimates for mu
% and sigma
function [X,cor,mu_bi,sigma_bi,mu_hat,sigma_hat,x]=SingleFileAnalysis(fid,mu_guess,sigma_guess,subID,TestType)
% 
% test=3;
% if test == 1
%    teststr=' Visual';
% elseif test == 2
%    teststr=' Auditory';
% elseif test ==3
%    teststr=' Tactlie';
% else
%    teststr=' Vestibular';
% end

%% Extract Data
guess_rate = 0.5;

master = importdata(fid);
if TestType == 2
X = master.data.Sheet1(:,2);
cor = master.data.Sheet1(:,4);
else
X = master.data.Sheet1(:,2);
X = abs(X);
cor = master.data.Sheet1(:,4);
end

% Fit binary threshold
x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mu_guess, sigma_guess]);
mu_bi = x(1);
sigma_bi = x(2);

% Fit lapse threshold
lambda_max = 0.06;
x = fminsearch(@(x) two_int_fit_lapse(x, X, cor, lambda_max), [0 mu_guess sigma_guess]);
lambda_hat = x(1);
mu_hat = x(2);
sigma_hat = x(3);

% %% Create 3 Panel Plot
% titlestr=[num2str(subID) teststr];
% xstr='Stimulus'; ystr='Percent Correct';
% if length(X) == 100
%     clusterSize = 20;
% elseif length(X) == 50
%     clusterSize = 10;
% end
% [x,prob_bi,prob_hat]=plotSBC(X,cor,lambda_hat,mu_bi,sigma_bi,mu_hat,sigma_hat,guess_rate,titlestr,xstr,ystr,clusterSize);