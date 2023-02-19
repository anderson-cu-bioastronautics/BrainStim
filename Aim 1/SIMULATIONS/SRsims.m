

clear all; close all; clc; 

tic

mu_mean = 5;
mu_sd = 1;
sig_mean = 1;
sig_sd = 0.2;
mu = normrnd(mu_mean, mu_sd);
sigma = normrnd(sig_mean, sig_sd);       

plot_count=5;
for M=1:plot_count   
%% Set up loop for different trial numbers
K = 6; %Number of data points per subject at each noise level 
mu_bi_vect=zeros(1,K);
for n=1:K 
noiselevels = [0 10 20 30 40 50];    
L=noiselevels(n);
            %% Parameters of the "theoretical" subject's underlying values
            %sigma = 1;  % sigma value, "slope" of psychometric curve
            %mu = 5;   % underlying "threshold" at which subject gets 75% correct
            guess_rate = 0.5;  % the level at which subject guesses correctly, normally 0.5
            initial_stim=10;
            
            
            % calculate the stimuli in which the subject would get 79.4% correct
            % theoretically, given the underlying mu and sigma
            thresh_794 = icdf('norm', (0.794 - guess_rate)/(1-guess_rate), mu, sigma);
            %% Perform the 3 down 1 up staircase, simulating the subject's responses
            N = 30;       % number of trials in the simulated session
            
            lapse_rate = 0; % we assuming the underlying lapse rate is 0 (i.e., the subject never lapses)
            %initial_stim = 10; %staircase
            plot_on = 0;
            plot_on2 = 0;
            % pest_mod_2int performs the 3 down 1 up staircase for a 2 interval
            % detection task
            [X, Y, cor, lapses] = pest_mod_2int(N, mu, sigma, guess_rate, lapse_rate, plot_on, initial_stim);
            % lapses = 1 if lapses, 0 if not. since we are assuming no lapses, these
            % will all be 0's
            X = X';             % the simuli levels, negative values correspond to first interval and positive to second inteval, the absolute value is the magnitude of the stimulus
            Ybi = [~Y'; Y'];    % the subject's responses. 0 = respond stim is in first interval, 1 = respond stim is in second interval
            cor = cor';         % cor = 0 if incorrect, = 1 if the responses is correct
            
            %% Fit the binary psychometric curve using the simple Gaussian from 0.5 to 1
            mu_guess = 4;
            sigma_guess = 1;
            x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mu_guess, sigma_guess]);
            mu_bi = x(1);
            sigma_bi = x(2);
            thresh_794_bi = icdf('norm', (0.794 - guess_rate)/(1-guess_rate), mu_bi, sigma_bi);
            
            
            %% Fit the binary psychometric curve with a guess rate differing from 0.5 and a non-zero lapse rate
            lambda_max = 0.06;
            gamma_max = 0.06;
            x = fminsearch(@(x) two_int_fit_guesslapse(x, X, cor, lambda_max, gamma_max), [0 0 mu_guess sigma_guess]);
            lambda_hat = x(1);
            gamma_hat = x(2);
            mu_hat = x(3);
            sigma_hat = x(4);
            thresh_794_hat = icdf('norm', (0.794 - guess_rate - gamma_hat)/(1-0.5-gamma_hat-lambda_hat), mu_hat, sigma_hat);
            
            %% Plot the results
            % Set up a x-vector for fitting
            xmax = max(abs(X));
            xval = linspace(0, xmax);
            
            % Cumulative normal distribution of the fitted psychometric curve
            prob_bi = guess_rate + (1-guess_rate-0)*cdf('norm', xval, mu_bi, sigma_bi);
            
            % Cumulative normal distribution of the fitted psychometric curve with
            % a guess rate differing from 0.5 and non-zero lapse rate
            prob_hat = guess_rate + gamma_hat + (1-0.5-gamma_hat-lambda_hat)*cdf('norm', xval, mu_hat, sigma_hat);
            
            % Cumulative normal distribution of the underlying psychometric curve
            % (note: this is never known in a real subject)
            prob = guess_rate + (1-guess_rate-0)*cdf('norm', xval, mu, sigma);
            
            if plot_on2==1
                figure;
                plot(xval, prob, '--'); hold on; % plot the underlying psychometric curve
                plot(xval, prob_bi);            % plot the fitted psychometric curve (simple)
                plot(xval, prob_hat, '--');   % plot the fitted psychometric curve (guess and lapse rates)
            else
            end
            
            % plot the mean probabilities with the size of the marker representing the
            % number of trials at that particular stimulus level
            [Xunique, ix, ixu] = unique(X); % find the unique stimuli levels tested
            % Average the Y outcomes at each unique X value
            sumUnique = zeros(2,length(Xunique));
            Punique = zeros(2,length(Xunique));
            Lunique = zeros(1,length(Xunique));
            for i = 1:length(Xunique)
                coratXunique = cor(:,ixu == i);          % find the correct outcomes for the ith unique X value
                sumUnique(:,i) = sum(coratXunique,2);   % total number of correct responses for the ith unique X value
                Lunique(i) = size(coratXunique,2);    % find the number of trials for the ith unique X value
                Punique(:,i) = sumUnique(:,i)/Lunique(i);      % find the probabilities for a correct response at the ith unique X value
            end
            for i = 1:length(Xunique)
                % find the marker size with a max of 12 and a min of 4
                if max(Lunique) == min(Lunique)
                    msize = 8;
                else
                    msize = (12-4)/(max(Lunique)-min(Lunique)) * (Lunique(i)-min(Lunique)) + 4;
                end
                
                if plot_on2==1
                    plot(Xunique(i), Punique(2,i), 'ko', 'MarkerSize', msize); % plot marker with the correct size
                else
                end
                %     if SPMdata(2,i)>= 5 && SPMdata(3,i) >=5  % need at least 5 of each type of response to use this approximation of CI, otherwise have to use exact methods, which I have not implemented yet
                %         plot(Xunique(i)*[1 1], Punique(2,i)+1.96*sqrt(Punique(2,i)*(1-Punique(2,i))/Lunique(i))*[-1 1], 'k');  % plot 95% confidence intervals
                %     end
            end
            
            if plot_on2==1
                legend('underlying', 'fitted simple', 'fitted guess/lapse', 'data', 'Location', 'SouthEast')
                xlabel('Stimulus (e.g., tactile vibration magnitude)'); ylabel('Likelihood of Correct Response')
                ylim([0 1.05])
            else
            end
mu_bi_vect(1,n)=mu_bi;  
end
%%
subplot(3,2,M)
plot(noiselevels,mu,'r-*',noiselevels,mu_bi_vect,'-o')
grid on
title(['Simulated Subject, No SR, N=' num2str(N)])
xlabel('Noise Level')
ylabel('Threshold') 
ylim([0 10])
end
toc