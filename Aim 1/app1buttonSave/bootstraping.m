function [musi_bo] = bootstraping(mu_bi,sigma_bi,mu_hat,sigma_hat,mu_guess,sigma_guess,TestType)

%% Perform the 3 down 1 up staircase, simulating the subject's responses
sim=200;
for i=1:sim
    guess_rate = 0.5;  % the level at which subject guesses correctly, normally 0.5
    lapse_rate = 0; % we assuming the underlying lapse rate is 0 (i.e., the subject never lapses)
    plot_on = 0;
    %% Fit the binary psychometric curve
    if TestType == 1 % Visual simulation
    N = 50;       % number of trials in the simulated session
    initial_stim=0.5;
    % pest_mod_2int performs the 3 down 1 up staircase for a 2 interval
    % detection task
    [X, Y, cor, lapses] = pest_mod_2int(N, mu_bi, sigma_bi, guess_rate, lapse_rate, plot_on, initial_stim);
    elseif TestType == 2 % Auditory simulation
    N = 100;       % number of trials in the simulated session
    initial_stim=20;
    % pest_mod_2int performs the 3 down 1 up staircase for a 2 interval
    % detection task
    [X, Y, cor, lapses] = pest_mod_2int_Audio(N, mu_bi, sigma_bi, guess_rate, lapse_rate, plot_on, initial_stim);    
    elseif TestType == 3 % Tactile simulation
    N = 50;       % number of trials in the simulated session
    initial_stim=.4;
    % pest_mod_2int performs the 3 down 1 up staircase for a 2 interval
    % detection task
    [X, Y, cor, lapses] = pest_mod_2int(N, mu_bi, sigma_bi, guess_rate, lapse_rate, plot_on, initial_stim);    
        
    elseif TestType == 4 % Vestibular simulation
        
    end
    
    % lapses = 1 if lapses, 0 if not. since we are assuming no lapses, these
    % will all be 0's
    X = X';             % the simuli levels, negative values correspond to first interval and positive to second inteval, the absolute value is the magnitude of the stimulus
    Ybi = [~Y'; Y'];    % the subject's responses. 0 = respond stim is in first interval, 1 = respond stim is in second interval
    cor = cor';         % cor = 0 if incorrect, = 1 if the responses is correct
    
    x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mu_guess, sigma_guess]);
    mus_bi = x(1);
    sigmas_bi = x(2);
    thresh_794_bi = icdf('norm', (0.794 - guess_rate)/(1-guess_rate), mus_bi, sigmas_bi);
    
    %% Fit the lapse psychometric curve
    if TestType == 1 % Visual simulation
    N = 50;       % number of trials in the simulated session
    initial_stim=0.5;
    % pest_mod_2int performs the 3 down 1 up staircase for a 2 interval
    % detection task
    [X, Y, cor, lapses] = pest_mod_2int(N, mu_hat, sigma_hat, guess_rate, lapse_rate, plot_on, initial_stim);
    elseif TestType == 2 % Auditory simulation
    N = 100;       % number of trials in the simulated session
    initial_stim=20;
    % pest_mod_2int performs the 3 down 1 up staircase for a 2 interval
    % detection task
    [X, Y, cor, lapses] = pest_mod_2int_Audio(N, mu_hat, sigma_hat, guess_rate, lapse_rate, plot_on, initial_stim); 
    
    elseif TestType == 3 % Tactile simulation
        
    elseif TestType == 4 % Vestibular simulation
        
    end
    
%     if i == 1
%         plot(1:length(X),X,'r.','MarkerSize',20);
%     end
    lambda_max = 0.06;
%     gamma_max = 0.06;
    x = fminsearch(@(x) two_int_fit_lapse(x, X, cor, lambda_max), [0 mu_guess sigma_guess]);
    lambda_hat = x(1);
    mus_hat = x(2);
    sigmas_hat = x(3);
%     thresh_794_hat_boot = icdf('norm', (0.794 - guess_rate)/(1-0.5-lambda_hat), mus_hat, sigmas_hat);
    musi_all(i,:)=[mus_bi sigmas_bi mus_hat sigmas_hat];
    
end

% Take 95% CI and organize data
mu_se_bi_list=sort(musi_all(:,1));
sigma_se_bi_list=sort(musi_all(:,2));
mu_se_hat_list=sort(musi_all(:,3));
sigma_se_hat_list=sort(musi_all(:,4));
mu_se_bi=[mu_se_bi_list(ceil(sim*.025)) mu_se_bi_list(floor(sim*.975))];
sigma_se_bi=[sigma_se_bi_list(ceil(sim*.025)) sigma_se_bi_list(floor(sim*.975))];
mu_se_hat=[mu_se_hat_list(ceil(sim*.025)) mu_se_hat_list(floor(sim*.975))];
sigma_se_hat=[sigma_se_hat_list(ceil(sim*.025)) sigma_se_hat_list(floor(sim*.975))];

musi_bo=[mu_se_bi sigma_se_bi mu_se_hat sigma_se_hat];
