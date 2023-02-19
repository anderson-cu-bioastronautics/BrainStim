%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to observe collective data and analyze the personal statistics
% Started by Sage and Abby
%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Get the deets
Subject = input('What is the subject number : ');
Retest = input('Enter \n 0 for Original Test Analysis \n 1 for Retest Analysis \nInput: ');
currentFolder = pwd;

% Subject file information
if Subject >9
    Subject_str = num2str(Subject);
else
    Subject_str = ['0',num2str(Subject)];
end
subFolder = ['Subject_Data//Subject_',Subject_str];

% figure out testing order
Orderfile = ['TestingOrder_' num2str(Subject) '.xls'];
OrderFile = fullfile(currentFolder,'Subject_Data',['Subject_' num2str(Subject)],Orderfile);
LevelData = readtable(OrderFile);
SRlevels = table2array(LevelData); % add more later, probably use a test randomizer excel sheet like other aims


x4mult = 1; row_ind = 1;
for trial = 1:length(SRlevels)*4
    % load in file data
    if trial >9
        trial_str = num2str(trial);
    else
        trial_str = ['0',num2str(trial)];
    end
    file_str = ['Subject_' Subject_str '_Set0' num2str(x4mult) '_T0' num2str(row_ind) '_RAW.mat'];
    filename = fullfile([pwd '//' subFolder '//' file_str]);
  % [rms_total, fracTimeIn, audStop, tacStop, Choice_P, CAL_Score, totalHeadingChanges] = Post_trial_script_fun(filename, trial);
    [rms_total, fracTimeIn, audStop, Choice_P, CAL_Score, totalHeadingChanges] = Post_trial_script_fun(filename, trial);
    
    % Sort file information
    rms_total_vect(row_ind,x4mult)=rms_total;
    fracTimeIn_vect(row_ind,x4mult)=fracTimeIn;
    audStop_vect(row_ind,x4mult)=audStop(1);
    audStop_vect(row_ind+4,x4mult)=audStop(2);
 %   tacStop_vect(row_ind,x4mult)=tacStop(1);
  %  tacStop_vect(row_ind+4,x4mult)=tacStop(2);
    Choice_P_vect(row_ind,x4mult)=Choice_P;
    CAL_Score_vect(row_ind,x4mult)=CAL_Score;
    Head_Change_vect(row_ind,x4mult)=totalHeadingChanges;
    row_ind = row_ind+1;
     
    if floor(trial/(4*x4mult)) == 1 % help divide into SR columns
        x4mult = x4mult+1;
        row_ind = 1;
    else
    end
end

%% Ok so like... create some distributions and run with it, god I hope it's normal
% run some normality tests when you have an actual dataset

% Flight performance : Distance to reticle
rms_total_mu = mean2(rms_total_vect);
rms_total_std = std2(rms_total_vect);

% Flight performance : Time Input
fracTimeIn_mu = mean2(fracTimeIn_vect);
fracTimeIn_std = std2(fracTimeIn_vect);

% Perception Performance : Auditory
audStop_mu = mean2(audStop_vect);
audStop_std = std2(audStop_vect);

% % Perception Performance : Tactile
% tacStop_mu = mean2(tacStop_vect);
% tacStop_std = std2(tacStop_vect);

% Decision Performance : Ranked landings
Choice_P_mu = mean2(Choice_P_vect);
Choice_P_std = std2(Choice_P_vect);

% Decision Performance : Crash, Abort, Land
CAL_Score_mu = mean2(CAL_Score_vect);
CAL_Score_std = std2(CAL_Score_vect);

% Decision Performance : Crash, Abort, Land
Head_Change_mu = mean2(Head_Change_vect);
Head_Change_std = std2(Head_Change_vect);

%% Calculate the Probabilities of each SR level across the distribution
for i=1:length(SRlevels)
    % Component performance analysis
   rms_total_p(i) = normcdf(mean(rms_total_vect(:,i)),rms_total_mu,rms_total_std);
   fracTimeIn_p(i) = normcdf(mean(fracTimeIn_vect(:,i)),fracTimeIn_mu,fracTimeIn_std);
   audStop_p(i) = normcdf(mean(audStop_vect(:,i)),audStop_mu,audStop_std);
  % tacStop_p(i) = normcdf(mean(tacStop_vect(:,i)),tacStop_mu,tacStop_std);
   Choice_P_p(i) = normcdf(mean(Choice_P_vect(:,i)),Choice_P_mu,Choice_P_std);
   CAL_Score_p(i) = normcdf(mean(CAL_Score_vect(:,i)),CAL_Score_mu,CAL_Score_std);
   Head_Change_Score_p(i) = normcdf(mean(Head_Change_vect(:,i)),Head_Change_mu,Head_Change_std);
   
   %Total_Performance(i) = (1-rms_total_p(i)) + (1-fracTimeIn_p(i)) + (1-audStop_p(i)) + (1-tacStop_p(i)) + Choice_P_p(i) + CAL_Score_p(i) + (1-Head_Change_Score_p(i));
   Total_Performance(i) = (1-rms_total_p(i)) + (1-fracTimeIn_p(i)) + (1-audStop_p(i)) + Choice_P_p(i) + CAL_Score_p(i) + (1-Head_Change_Score_p(i));
end

%% Plot Data
size = length(SRlevels)
 XtickMarks = 1:size;
for i = 1:size
    if SRlevels(i) == 0
        Xlabel{i} = 0;
    elseif SRlevels(i) == 1
        Xlabel{i} = 'MM';
    else
        Xlabel{i} = SRlevels(i);
    end
end

plot(XtickMarks,Total_Performance); 
ylabel('Total Performance Score')
xlabel('SR Level') % Change this to condition later
xticks(XtickMarks); xticklabels(Xlabel); 

if ~Retest 

  for i = 1:length(SRlevels)  %Sort into VSR and ASR
        if SRlevels(i) == 0
            vsr(i) = 2;      
        elseif SRlevels(i) <= 1 %vsr
            vsr(i) = 1;
        else
            vsr(i) = 0;
        end
  end
    
    %Calculate and output VSR and ASR level with best feedback
    vsrFeed = Total_Performance(vsr == 1); %find max of feedback when vsr
    maxVSR = Total_Performance(max(vsrFeed) == Total_Performance); %find max of vsr 
    maxVSRLev = SRlevels(maxVSR == Total_Performance);

    asrFeed = Total_Performance(vsr == 0);
    maxASR = Total_Performance(max(asrFeed) == Total_Performance);
    maxASRLev = SRlevels(maxASR == Total_Performance);
    

    %ramdomize retest levels
    retestLevels = {'sham', 'MMSR', 'ASR', 'VSR'};
    randOrder=randperm(4);
    retestOrder = retestLevels(randOrder);

    %Print retest levels and order
    fprintf('\n')
    MMSRLevOut = ['Retest MMSR:   VSR = ', num2str(maxVSRLev), ' mA ', '   ASR = ', num2str(maxASRLev), ' dB'];
    disp(MMSRLevOut);
    retestOrderOut = ['Retest Order: ', char(retestOrder(1)), ', ', char(retestOrder(2)), ', ', char(retestOrder(3)), ', ', char(retestOrder(4))];
    disp(retestOrderOut);
    
    %
    %Write Retest Levels To file
    retestWrite = [0, 1, maxASRLev, maxVSRLev];
    retestForFile = retestWrite(randOrder);
    xlStart = 'H1'; %first column of data to write to assuming 7 original tests 
    xlswrite(OrderFile, retestForFile, 1, xlStart)
    
   
end





