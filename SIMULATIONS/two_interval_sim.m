% This code simulates binary threshold data for a two-interval task (e.g.,
% in tactile thresholds this might be a detection task of whether the
% stimulus was in the first interval or the second interval. Similar for
% auditory). The data is simulated using a 3 down 1 up staircase. It is
% then fit with a cumulative gaussian psychometric curve that has been
% scaled to go from 0.5 to 1, yielding a mu and sigma. The mu parameter is
% the stimulus level at which the subject is getting 75% correct. We also
% calculate the level at which the subject is getting 79.4% correct
% (corresponding to the 3D1U staircase reversal level). One of these parameters
% (and not sigma) is typically defined as the "threshold"

clear all; close all; clc; 

tic

%% Set up loop for different sigma/staircase values
% for s=1:3
%     Svect=[5 10 20]; %staircase values 
%     Svect=[0.2 1 2.5]; %sigma values
%     sigma=Svect(s);
     sigma=1;
%     initial_stim = Svect(s); %staircase    
     initial_stim=10;
    
    %% Set up loop for different trial numbers
    K = 10; %number of different numbers of trials in Nvect
    stdevvect=zeros(1,K); %set up matrices to store K # of values in
    cvvect=zeros(1,K);
    cvvect2=zeros(1,K);
    for n=1:K
        Nvect=[15 20 25 30 35 40 50 75 100 125];
        N=Nvect(n); %trial number        
        
        %% Set up loop for number of simulations to run        
        M_total=1000; %number of simulations to run, set as even multiple of 100 to avoid line 147 error
        
        Mplot_on=0; %histogram plot for each set of M_total simulations
        
        sim_vect=zeros(1,M_total);
        sim_vect2=zeros(1,M_total);
        
        for M=1:M_total
            %% Parameters of the "theoretical" subject's underlying values
            %sigma = 1;  % sigma value, "slope" of psychometric curve
            mu = 5;   % underlying "threshold" at which subject gets 75% correct
            guess_rate = 0.5;  % the level at which subject guesses correctly, normally 0.5
            
            % calculate the stimuli in which the subject would get 79.4% correct
            % theoretically, given the underlying mu and sigma
            thresh_794 = icdf('norm', (0.794 - guess_rate)/(1-guess_rate), mu, sigma);
            %% Perform the 3 down 1 up staircase, simulating the subject's responses
            %N = 50;       % number of trials in the simulated session
            
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

% Add after finishing staircase and sigma sims to replace above section
% Download new 2_int_fit from email
%             lambda_max = 0.06;
%             gamma_hat = 0;
%             x = fminsearch(@(x) two_int_fit_lapse(x, X, cor, lambda_max), [0 mu_guess sigma_guess]);
%             lambda_hat = x(1);
%             mu_hat = x(2);
%             sigma_hat = x(3);
%             thresh_794_hat = icdf('norm', (0.794 - guess_rate - 0)/(1-0.5-0-lambda_hat), mu_hat,   sigma_hat);

            
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
            
            sim_vect(1,M)=mu_bi; %add mu_bi/mu_hat/thresh etc to the vector to plot later
            sim_vect2(1,M)=thresh_794_bi;
            M=M+1
        end
        
        if Mplot_on==1
            figure
            hist(sim_vect,M_total) %histogram of mu values
            xlabel('Estimated mu value')
            ylabel('Frequency')
        else
        end
        
        sorted_sim_vect = sort(sim_vect);     % sort the value from low to high
        lower_per = sorted_sim_vect(round(0.17*M_total));    % grab the 17th percentile one, for 10,000 sims, this is the 1,700 value starting from the lowest
        upper_per = sorted_sim_vect(round(0.83*M_total));   % grab the 83rd percentile one
        range_per = upper_per - lower_per;   % this is roughly analogous to the standard deviation
        CV_range_per = range_per / (2*mu);    % just like when calculating CV from stdev, normalize by the underlying value.
        
        sorted_sim_vect2 = sort(sim_vect2);     % sort the value from low to high
        lower_per2 = sorted_sim_vect2(round(0.17*M_total));    % grab the 17th percentile one, for 10,000 sims, this is the 1,700 value starting from the lowest
        upper_per2 = sorted_sim_vect2(round(0.83*M_total));   % grab the 83rd percentile one
        range_per2 = upper_per2 - lower_per2;   % this is roughly analogous to the standard deviation
        CV_range_per2 = range_per2 / (2*thresh_794);        
%              avgmu=mean(sim_vect); %avg of mu_bi values
%              stdevmu=std(sim_vect); %std dev of mu_bi values
%              cvmu=stdevmu/mu; %coefficient of variance
%         
%              avgmu2=mean(sim_vect2);
%              stdevmu2=std(sim_vect2);
%              cvmu2=stdevmu2/thresh_794;        
               
        cvvect(1,n)=CV_range_per;
        cvvect2(1,n)=CV_range_per2;
%              cvvect(1,n)=cvmu; %add cv values to an array for each N
%              cvvect2(1,n)=cvmu2; 
%     
    disp(N)
%    end
%%     
plot(Nvect,cvvect,'r-o',Nvect,cvvect2,'b-o')      
%     if s==1
%         line1=cvvect2;        
%     elseif s==2
%         line2=cvvect2;
%     else
%         line3=cvvect2;
%     end
    
 end
%%
%plot(Nvect,line1,'b-o',Nvect,line2,'b-*',Nvect,line3,'b-+')
grid on
title('Comparison of Mu-bi and Thresh-794-bi')
xlabel('# of Trials')
ylabel('Non-Parametric CV') %Non-parametric cv if using sort method
%legend('Initial stim=5','Initial stim=10','Initial stim=20','Location', 'NorthEast')
legend('Mu-bi','Thresh-bi-794','Location','NorthEast')
toc