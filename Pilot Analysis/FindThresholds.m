% Function to calculate and find thresholds for all conditions

function [AllData]=FindThresholds(subID,TestType,SRType,levels)

%% define files to be read and define initial threshold guesses
currentFolder=pwd;
typeFolder='RA1Data';
subFold=num2str(subID);
switch TestType
    case 1 % define visual folder information
        mu_guess = 0.2; sigma_guess = 0.05; N=50;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Visual_0mA_'];
            parseFold2=['dB.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Visual_'];
            parseFold2=['mA_0dB.xls'];
        else
        end
    case 2 % define auditory folder information
        mu_guess = 8; sigma_guess = 3; N=100;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Auditory_0mA_'];
            parseFold2=['dB.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Auditory_'];
            parseFold2=['mA_0dB.xls'];
        else
        end
    case 3 % define tactile folder information
        mu_guess = 0.3; sigma_guess = 0.1; N=50;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Tactile_0mA_'];
            parseFold2=['dB.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Tactile_'];
            parseFold2=['mA_0dB.xls'];
        else
        end
    case 4
        mu_guess = 0; sigma_guess = 0; N=100;
        if SRType == 1
            parseFold1=['BSRA1_ASR_' num2str(subID) '_Vestibular_0mA_'];
            parseFold2=['dB.xls'];
        elseif SRType == 2
            parseFold1=['BSRA1_VSR_' num2str(subID) '_Vestibular_'];
            parseFold2=['mA_0dB.xls'];
        else
        end
end

count=1;
%% Loop through levels, calculate thresholds, simulate oh my!
guess_rate = 0.5;  % the level at which subject guesses correctly, normally 0.5
for i=1:length(levels)
    
    % Determine the file you are entering... slightly ugly here
    level=levels(i);
    if sign(level) == -1 && floor(level) == level % folder definitions for negative SR levels
        fileID=[parseFold1 'neg' num2str(abs(level)) parseFold2];
        fid=fullfile(currentFolder,typeFolder,subFold,fileID);
        master=importdata(fid);
    elseif sign(level) == 1 && floor(level) ~= level % definitions for non-integer SR levels
        fileID=[parseFold1 num2str(floor(level)) 'p' num2str(10*rem(level,1)) parseFold2];
        if TestType == 2
            fileID=[parseFold1 num2str(floor(level)) 'p' num2str(10*rem(level,1)) parseFold2];
        end
        fid=fullfile(currentFolder,typeFolder,subFold,fileID);
        master=importdata(fid);
    elseif sign(level) == -1 && floor(level) ~= level % folder definitions for negative SR levels
        fileID=[parseFold1 'neg' num2str(abs(ceil(level))) 'p' num2str(abs(10*rem(level,1))) parseFold2];
        fid=fullfile(currentFolder,typeFolder,subFold,fileID);
        master=importdata(fid);
    else % definition for all whole number SR levels
        fileID=[parseFold1 num2str(abs(level)) parseFold2];
        fid=fullfile(currentFolder,typeFolder,subFold,fileID)
        master=importdata(fid);
    end
    
    % separate out the data, do NOT take absolute value for auditory
    % this was already accounted for in the test code for auditory
    if TestType == 2
        X = master.data.Sheet1(:,2);
        cor = master.data.Sheet1(:,4);
    else
        X = master.data.Sheet1(:,2);
        X = abs(X);
        cor = master.data.Sheet1(:,4);
    end
    
    % Fit binary thresholds and find their jack knife bars
    disp(['Data binary, noise is ' num2str(level)])
    x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mu_guess, sigma_guess]);
    mu_bi = x(1);
    sigma_bi = x(2);
    %     thresh_794_bi = icdf('norm', (0.794 - guess_rate)/(1-guess_rate), mu_bi, sigma_bi);
    disp(['Jackknife binary, noise is ' num2str(level)])
    [musi_bi] = jackknife(@just794_bi,X,cor,mu_guess,sigma_guess);
    mu_se_bi = sqrt( (N-1)/N*sum( (mu_bi - musi_bi(:,1)).^2) );
    sigma_se_bi = sqrt( (N-1)/N*sum( (sigma_bi - musi_bi(:,2)).^2 ) );
    
    % Fit lapse thresholds and find their jack knife bars
    lambda_max = 0.06;
    %     gamma_max = 0.06;
    disp(['Data lapse, noise is ' num2str(level)])
    x = fminsearch(@(x) two_int_fit_lapse(x, X, cor, lambda_max), [0 mu_guess sigma_guess]);
    lambda_hat = x(1);
    mu_hat = x(2);
    sigma_hat = x(3);
    %     thresh_794_hat = icdf('norm', (0.794 - guess_rate)/(1-0.5-lambda_hat), mu_hat, sigma_hat);
    
    disp(['Jackknife lapse, noise is ' num2str(level)])
    [musi_hat] = jackknife(@just794_hat,X,cor,mu_guess,sigma_guess);
    mu_se_hat = sqrt( (N-1)/N*sum( (mu_hat - musi_hat(:,1)).^2) );
    sigma_se_hat = sqrt( (N-1)/N*sum( (sigma_hat - musi_hat(:,2)).^2 ) );
    
    disp(['Bootstrap both, noise is ' num2str(level)])
    rng default;
    % Now bootstrap our fine data here, 2000 simulations (shorten it in function if necessary)
    [musi_bo] = bootstraping(mu_bi,sigma_bi,mu_hat,sigma_hat,mu_guess,sigma_guess,TestType);
    
    % shove and save our data in cell arrays
    AllData{i,1} = X; % stim levels
    AllData{i,2} = cor; % correct?
    AllData{i,3} = [mu_bi sigma_bi]; % mus
    AllData{i,4} = [mu_hat sigma_hat]; % sigmas
    AllData{i,5} = [mu_se_bi sigma_se_bi]; % binary jack knife
    AllData{i,6} = [mu_se_hat sigma_se_hat]; % lapse jack knife
    AllData{i,7} = [musi_bo(1:2) musi_bo(3:4)]; % binary bootstrap
    AllData{i,8} = [musi_bo(5:6) musi_bo(7:8)]; % lapse bootstrap
    
    % Grab probability data for Psych Curves
    xmax = max(abs(X));
    xval = linspace(0, xmax);
    
    % Cumulative normal distribution of the fitted psychometric curve
    prob_bi = guess_rate + (1-guess_rate-0)*cdf('norm', xval, mu_bi, sigma_bi);
    
    % Cumulative normal distribution of the fitted psychometric curve with
    % a guess rate differing from 0.5 and non-zero lapse rate
    prob_hat = guess_rate + (1-0.5-lambda_hat)*cdf('norm', xval, mu_hat, sigma_hat);
    
    AllData{i,9}=xval;
    AllData{i,10}=prob_bi;
    AllData{i,11}=prob_hat;
    AllData{i,12}=lambda_hat;
    
end