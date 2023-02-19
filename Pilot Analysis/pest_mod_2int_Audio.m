% clear all;clc;close all

function [tstim, firstorsecond, cor, lapses] = pest_mod_2int_Audio(trial_no, mu_air, sigma_air, guess_rate, lapse_rate, plot_on, initial_stim)

initial_stim=10^initial_stim;

sim_no=1; 

n_down_PEST=3;
n_up_PEST=1;

%% parameter initialization

%min_delta_log=log10(2^(1/32));
min_delta_log=.5;
max_delta_log=10;

correct=[]; first_or_second=[]; trial_stim=[];

%% simulation main part
for sim_cnt=1:sim_no

    k=0;   r_cnt=0;
    n_down_init_cnt=0;  n_up_init_cnt=0;   
    n_down_PEST_cnt=0; 	n_up_PEST_cnt=0;
    double_down=0;      double_up=0;    %Initially no doubling down
    steps=1;                            %Upon first entering PEST trials, stim level will have stepped up one level
    
	stim_factor=max_delta_log;
	delta_log_stim=((stim_factor));

    %initialize variables
    n_down_init_cnt=0;  n_up_init_cnt=0;  
    n=0;  %Initialize n, where n is the trial number (1, 2, 3, .... N)
    lapses = zeros(trial_no, 1);

    %Initially set errors equal to zero
    reversal=0;
    stim_level=initial_stim;    
    stim_level_log=log10(stim_level);                
    up_or_down=-1;         

    while (k < trial_no )   %This loop is for the initial stimuli until the first mistake
        j=0; order=[];
        while (n_down_PEST_cnt < n_down_PEST && n_up_PEST_cnt < n_up_PEST )
            j=j+1;
            n=n+1;
            k=k+1;                   
            order(j)=(randi(2,1)*2-3);

%             trial_stim(n)= stim_level * order(j);  
            trial_stim(n) = stim_level_log;

%             p_system=lambda_init+(1-2*lambda_init)*0.5*(1 + erf((trial_stim(n)-mu_air)/(sqrt(2)*sigma_air)));  %This line with following if loop is a simulation of the actual experiment for + stimuli
            p_system = guess_rate+(1-guess_rate-lapse_rate)*0.5*(1 + erf((trial_stim(n)-mu_air)/(sqrt(2)*sigma_air)));  %This line with following if loop is a simulation of the actual experiment
            
%             lapses(n) = random('binomial', 1, lapse_rate);
            lapses(n) = binornd(1,lapse_rate);
            if lapses(n)
                p_system = 0.5;
            end

            
%             if sign(trial_stim(n)*(rand(1)-p_system))<0
            if (rand(1) - p_system) < 0
%                 l_or_r(n)=0.5*order(j)+0.5; 
                first_or_second(n) = 0.5*order(j) + 0.5;
                correct(n)=1;

                n_down_PEST_cnt=n_down_PEST_cnt+1;
                n_up_PEST_cnt=0;

                if (n_down_PEST_cnt == n_down_PEST)
                    if (up_or_down == 1)
                        reversal=reversal+1;
                        r_cnt=r_cnt+1;
                        delta_log_stim=delta_log_stim/2; if(delta_log_stim<min_delta_log),delta_log_stim=min_delta_log;end
                        steps=1;
                        double_down=0;
                        up_or_down=-1;
                    elseif(steps>2)  %if 4th time or greater always double it
                        delta_log_stim=delta_log_stim*2; if(delta_log_stim>max_delta_log),delta_log_stim=max_delta_log;end
                        double_down=1;
                        steps=steps+1;
                    elseif(steps==2 && double_down==0)  %if 3rd time and didn't just double, then double it
                        delta_log_stim=delta_log_stim*2; if(delta_log_stim>max_delta_log),delta_log_stim=max_delta_log;end
                        double_down=1;
                        steps=steps+1;
                    else
                        steps=steps+1;
                    end                            
                end
            else
%                 l_or_r(n)=-0.5*order(j)+0.5;
                first_or_second(n) = -0.5*order(j)+0.5;
                correct(n)=0;
                n_down_PEST_cnt=0;
                n_up_PEST_cnt=n_up_PEST_cnt+1;
                if (n_up_PEST_cnt == n_up_PEST)
                    if (up_or_down == -1)
                        reversal=reversal+1;
                        r_cnt=r_cnt+1;
                        delta_log_stim=delta_log_stim/2; if(delta_log_stim<min_delta_log),delta_log_stim=min_delta_log;end
                        steps=1;
                        double_up=0;
                        up_or_down=1; 
                    elseif(steps>2)  %if 4th time or greater always double it
                        delta_log_stim=delta_log_stim*2; if(delta_log_stim>max_delta_log),delta_log_stim=max_delta_log;end
                        double_up=1;
                        steps=steps+1;
                    elseif(steps==2 && double_up==0)  %if 3rd time and didn't just double, then double it
                        delta_log_stim=delta_log_stim*2; if(delta_log_stim>max_delta_log),delta_log_stim=max_delta_log;end
                        double_up=1;
                        steps=steps+1;
                    else
                        steps=steps+1;
                    end
                end
            end
        end

        n_down_PEST_cnt=0;  n_up_PEST_cnt=0;   

        if(up_or_down==1)    %When mth mistake is made
            stim_level_log=stim_level_log+delta_log_stim;
        else
            stim_level_log=stim_level_log-delta_log_stim;
        end
        stim_level=10^stim_level_log;  
    end
    
    cor(:,sim_cnt)=correct(1:trial_no)';
    firstorsecond(:,sim_cnt)=first_or_second(1:trial_no)';
%     lor(:,sim_cnt)=l_or_r(1:trial_no)';
    tstim(:,sim_cnt)=trial_stim(1:trial_no)';
    lapses = lapses(1:trial_no);

%     [b] = brglmfit(trial_stim,l_or_r','binomial','link','probit');
%     fit_brglm(1,sim_cnt)=-b(1)/b(2);    % mu
%     fit_brglm(2,sim_cnt)=1/b(2);        % sigma
% 
%     [b] = glmfit(trial_stim,l_or_r','binomial','link','probit');
%     fit_glm(1,sim_cnt)=-b(1)/b(2);      % mu 
%     fit_glm(2,sim_cnt)=1/b(2);          % sigma
    
    % Plot the sequence of responses
    if plot_on
        trials = 1:trial_no;
        ind_corr_first = (correct(1:trial_no) == 1).*(first_or_second(1:trial_no) == 0);
        ind_corr_second = (correct(1:trial_no) == 1).*(first_or_second(1:trial_no) == 1);
        ind_incorr_first = (correct(1:trial_no) ~= 1).*(first_or_second(1:trial_no) == 0);
        ind_incorr_second = (correct(1:trial_no) ~= 1).*(first_or_second(1:trial_no) == 1);
        
        figure; hold on;
        plot(trials(ind_corr_second==1), trial_stim(ind_corr_second==1), 'ro', 'MarkerSize', 4)
        plot(trials(ind_corr_first==1), trial_stim(ind_corr_first==1), 'ko', 'MarkerSize', 4)
        plot(trials(ind_incorr_second==1), trial_stim(ind_incorr_second==1), 'rx', 'MarkerSize', 6)
        plot(trials(ind_incorr_first==1), trial_stim(ind_incorr_first==1), 'kx', 'MarkerSize', 6)
        %     plot(trials(excluded_trials), abs(trial_stim(excluded_trials)), 'rs', 'MarkerSize', 10)
        xlabel('Trial Number');
        ylabel('Simulus Intensity, dB');
        % ylim([0 max(abs(trial_stim))]);
        text(0.25*trial_no, 0.9*max(trial_stim), 'o = first interval, correct response')
        text(0.25*trial_no, 0.85*max(trial_stim), 'o = second interval, correct response', 'Color', 'r')
        text(0.25*trial_no, 0.8*max(trial_stim), 'x = first interval, incorrect response')
        text(0.25*trial_no, 0.75*max(trial_stim), 'x = second interval, incorrect response', 'Color', 'r')
        box on;
    end
    
%     figure;
    
    
    clear correct first_or_second trial_stim b
end

% toc
