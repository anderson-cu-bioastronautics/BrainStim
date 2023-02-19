
function [] = wrapper_code_Tactile(app)
%% parameter initialization

n_down_PEST = 3;                % PEST staircase set up
n_up_PEST=1;

mu_guess = 0.3; sigma_guess = 0.1;

initial_stim = 1;

use_mags = 1;   % flag to use fixed magnitudes, discrete intervals
plot_on = 0;    % do you want to plot at the end, 1=yes, 0=no
fit_on = 1;     % do you want to calculate the threshold at the end, 1=yes, 0=no\

trial_no = 50;

mags=linspace(0.004, 1, 254);

%min_delta_log=log10(2^(1/32));
min_delta_log=log10(2^(1/8));
max_delta_log=log10(2^2);

correct=[]; l_or_r=[]; trial_stim=[];

%% simulation main part
% for sim_cnt=1:sim_no

k=0;   r_cnt=0;
n_down_init_cnt=0;  n_up_init_cnt=0;
n_down_PEST_cnt=0; 	n_up_PEST_cnt=0;
double_down=0;      double_up=0;    %Initially no doubling down
steps=1;                            %Upon first entering PEST trials, stim level will have stepped up one level

stim_factor=2;
delta_log_stim=(log10(stim_factor));

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
        app.n=n;
        app.nLabel.Text=num2str(n);
        k=k+1;
        order(j)=(randi(2,1)*2-3);

        if stim_level > 1
            stim_level = 1;
        end
        
        trial_stim(n)= stim_level * order(j);
        
        if use_mags
            % find the closest available stimulus and use that
            [dif, ind] = min(abs(abs(trial_stim(n)) - mags));
            trial_stim_run(n) = mags(ind) * order(j);
        else
            % just use the raw stimulus value
            trial_stim_run(n) = trial_stim(n);
        end

        
        % Tell the operator what they should run
        if n > trial_no
            sca;
            close all
            break
        elseif n==1
            uiwait
            app.Lamp.Color=[1 0 0];
            app.StimDisplayLabel.Text=sprintf('%g',abs(trial_stim_run(n)));
            [response] = hardwarefun_Tactile(trial_stim_run(n),app);
            l_or_r(n) = response;
        else
            %not complete
            app.StimDisplayLabel.Text=sprintf('%g',abs(trial_stim_run(n)));
            [response] = hardwarefun_Tactile(trial_stim_run(n),app);
            l_or_r(n) = response;
        end
        
        
        % check if the response was correct or incorrect and calculate the
        % magnitude of the next trial using PEST rules
        
        
        if l_or_r(n)*trial_stim_run(n) == 0
            n=n-1;
            k=k-1;
        elseif l_or_r(n)*trial_stim_run(n) > 0
            % CORRECT!
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
            % INCORRECT!
            correct(n)=0;
            app.correct(n) = correct(n);
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

% Pull together the outputs
trials = (1:trial_no)';
cor = correct(1:trial_no)';
lor = (l_or_r(1:trial_no)==1)';
tstim = trial_stim(1:trial_no)';
tstim_run = trial_stim_run(1:trial_no)';

output = [trials tstim tstim_run lor cor];


%% Fit the data
mu_guess=.5; sigma_guess=.1;
if fit_on
    [musi_bi]=just794_bi(abs(trial_stim_run),correct,mu_guess,sigma_guess);
    mu_est = musi_bi(1);
    sigma_est = musi_bi(2);
    display(['Threshold = ', num2str(mu_est), ' percentage displacement sigma =  ', num2str(sigma_est) ])
end
%% Plot the data
% Plot the sequence of responses
if plot_on
    % Plot what was actually run
    ind_corr_pos = (cor == 1).*(sign(tstim_run) == 1);
    ind_corr_neg = (cor == 1).*(sign(tstim_run) == -1);
    ind_incorr_pos = (cor ~= 1).*(sign(tstim_run) == 1);
    ind_incorr_neg = (cor ~= 1).*(sign(tstim_run) == -1);
    
    ylimit = max([max(abs(tstim_run)) max(abs(tstim))]);
    
    figure;
    subplot(3,1,1);
    hold on;
    plot(trials(ind_corr_pos==1), abs(tstim_run(ind_corr_pos==1)), 'ro', 'MarkerSize', 4)
    plot(trials(ind_corr_neg==1), abs(tstim_run(ind_corr_neg==1)), 'ko', 'MarkerSize', 4)
    plot(trials(ind_incorr_pos==1), abs(tstim_run(ind_incorr_pos==1)), 'rx', 'MarkerSize', 6)
    plot(trials(ind_incorr_neg==1), abs(tstim_run(ind_incorr_neg==1)), 'kx', 'MarkerSize', 6)
    xlabel('Trial Number');
    ylabel('Simulus Magnitude (deg)');
    ylim([0 ylimit]);
    title('Stimuli Run')
    box on;
    
    % Plot what we wanted to run
    ind_corr_pos = (cor == 1).*(sign(tstim) == 1);
    ind_corr_neg = (cor == 1).*(sign(tstim) == -1);
    ind_incorr_pos = (cor ~= 1).*(sign(tstim) == 1);
    ind_incorr_neg = (cor ~= 1).*(sign(tstim) == -1);
    
    subplot(3,1,2);
    hold on;
    plot(trials(ind_corr_pos==1), abs(tstim(ind_corr_pos==1)), 'ro', 'MarkerSize', 4)
    plot(trials(ind_corr_neg==1), abs(tstim(ind_corr_neg==1)), 'ko', 'MarkerSize', 4)
    plot(trials(ind_incorr_pos==1), abs(tstim(ind_incorr_pos==1)), 'rx', 'MarkerSize', 6)
    plot(trials(ind_incorr_neg==1), abs(tstim(ind_incorr_neg==1)), 'kx', 'MarkerSize', 6)
    xlabel('Trial Number');
    ylabel('Simulus Magnitude (deg)');
    ylim([0 ylimit]);
    title('Stimuli Desired');
    box on;
    
    
    % Average the Y outcomes at each unique X value
    [Xunique, ix, ixu] = unique(tstim_run);
    Punique = zeros(1,length(Xunique));
    Lunique = zeros(1,length(Xunique));
    for k = 1:length(Xunique)
        YatXunique = lor(ixu == k); % find the Y outcomes for the jth unique X value
        Lunique(k) = length(YatXunique);    % find the number of trials for the ith unique X value
        Punique(k) = mean(YatXunique);  % find the probability at the ith unique X value
    end
    
    
    subplot(3,1,3);
    hold on;
    X_vect = linspace(-ylimit, ylimit);
    prob_hat = cdf('norm', X_vect, mu_est, sigma_est);
    plot(X_vect, prob_hat, 'Color', 0.5*ones(1,3), 'LineWidth', 2)
    
    % plot the mean probabilities with the size of the marker representing the
    % number of trials
    for k = 1:length(Xunique)
        % find the marker size with a max of 12 and a min of 4
        if max(Lunique) == min(Lunique)
            msize = 8;
        else
            msize = (12-4)/(max(Lunique)-min(Lunique)) * (Lunique(k)-min(Lunique)) + 4;
        end
        plot(Xunique(k), Punique(k), 'ko', 'MarkerSize', msize);
    end
    legend('est curve', 'subject responses', 'Location', 'SouthEast')
    xlabel('Stimulus (deg)'); ylabel('Likelihood of Rightward Response')
    box on;
    
  %  name = [ 'PEST Tactile Pilot Testing Subject ' num2str(subject_number) '.jpg' ];
    %saveas(figure, name);
end

end