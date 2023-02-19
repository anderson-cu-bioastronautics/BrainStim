%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to observe collective data and analyze the personal statistics
% Started by Sage and Abby
%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
%close all; clc;

%% Get the deets
Subject = input('What is the subject number : ');
Retest = 1;
currentFolder = pwd;

% Subject file information
if Subject >=10
    Subject_str = num2str(Subject);
else
    Subject_str = ['0',num2str(Subject)];
end
subFolder = ['Subject_Data//PerformanceLZ//Subject_',Subject_str];

% figure out testing order
Orderfile = ['TestingOrder_' Subject_str '.xls'];
OrderFile = fullfile(currentFolder,subFolder,['Subject_' Subject_str],Orderfile);
LevelData = xlsread(OrderFile);
SRlevels = (LevelData); % add more later, probably use a test randomizer excel sheet like other aims
LocateFile = fullfile(currentFolder,subFolder,['Subject_' Subject_str]);

if any(Retest) == 1
  %  SRlevels = SRlevels(7:end);
else
    SRlevels = SRlevels(1:6);
end

x4mult = 1; row_ind = 1;
nantrial = zeros(length(SRlevels(7:end))*4, 1);
for trial = (1:length(SRlevels(7:end))*4) + 18
    % load in file data
    if trial >9
        trial_str = num2str(trial);
    else
        trial_str = ['0',num2str(trial)];
    end
    if x4mult < 4
        file_str = ['Subject_' Subject_str '_Set0' num2str(x4mult+6) '_T0' num2str(row_ind) '_RAW.mat'];
    else
        file_str = ['Subject_' Subject_str '_Set' num2str(x4mult+6) '_T0' num2str(row_ind) '_RAW.mat'];
    end
    filename = fullfile([pwd '//' subFolder '//' file_str]);
    if isfile(filename)                                         % If file exists, run post-trial script
        [rms_total, fracTimeIn, audStop, tacStop, Choice_P, CAL_Score, totalHeadingChanges] = Post_trial_script_fun(filename, trial, Retest);
        %[rms_total, fracTimeIn, audStop, Choice_P, CAL_Score, totalHeadingChanges] = Post_trial_script_fun(filename, trial);
    else                                                        % Otherwise, if data is missing, assign NaNs to metrics
        rms_total = NaN;
        fracTimeIn = NaN;
        audStop = NaN(2, 1);
        tacStop = NaN(2, 1);
        Choice_P = NaN;
        CAL_Score = NaN;
        totalHeadingChanges = NaN;
        nantrial(trial - 18) = 1;                               % Keep track of trials with missing data, since FractionalRankings() will drop NaNs
    end
    
    % Sort file information
    rms_total_mat(row_ind,x4mult)=rms_total;
    fracTimeIn_mat(row_ind,x4mult)=fracTimeIn;
    audStop_mat(row_ind,x4mult)=audStop(1);
    audStop_mat(row_ind+4,x4mult)=audStop(2);
    tacStop_mat(row_ind,x4mult)=tacStop(1);
    tacStop_mat(row_ind+4,x4mult)=tacStop(2);
    Choice_P_mat(row_ind,x4mult)=Choice_P;
    CAL_Score_mat(row_ind,x4mult)=CAL_Score;
    Head_Change_mat(row_ind,x4mult)=totalHeadingChanges;
    row_ind = row_ind+1;
    
    if floor((trial+2)/(4*(x4mult+5))) == 1 % help divide into SR columns
        x4mult = x4mult+1;
        row_ind = 1;
    else
    end
end

% save all this data
    SRPerfIndData = {rms_total_mat, fracTimeIn_mat, Head_Change_mat, audStop_mat, tacStop_mat, Choice_P_mat, CAL_Score_mat};
    savefilename = ['IndividualPerformanceData' Subject_str 'Retest.mat'];
    saveFilename = fullfile(LocateFile,savefilename);
    save(saveFilename,'SRPerfIndData');
    
%% Ok so like... rank them all in a matrix
% Sort ranking information
rms_array = reshape(rms_total_mat,[numel(rms_total_mat),1]);
fracTimeIn_array = reshape(fracTimeIn_mat,[numel(fracTimeIn_mat),1]);
audStop_array = reshape(audStop_mat,[numel(audStop_mat),1]);
tacStop_array = reshape(tacStop_mat,[numel(tacStop_mat),1]);
Choice_P_array = reshape(Choice_P_mat,[numel(Choice_P_mat),1]);
CAL_Score_array = reshape(CAL_Score_mat,[numel(CAL_Score_mat),1]);
Head_Change_array = reshape(Head_Change_mat,[numel(Head_Change_mat),1]);
inantrial = find(nantrial);                                 % Find indices of trials with missing data

% Compute fractional rankings
[rms_ranking] = FractionalRankings(rms_array);
[TimeIn_ranking] = FractionalRankings(fracTimeIn_array);
[audStop_ranking] = FractionalRankings(audStop_array);
[tacStop_ranking] = FractionalRankings(tacStop_array);
[Choice_P_ranking] = FractionalRankings(Choice_P_array);
[CAL_ranking] = FractionalRankings(CAL_Score_array);
[Head_Change_ranking] = FractionalRankings(Head_Change_array);

% Reinsert NaNs
for i = 1:length(inantrial)
    rms_ranking = [rms_ranking(1:inantrial(i) - 1); NaN; rms_ranking(inantrial(i):end)];
    TimeIn_ranking = [TimeIn_ranking(1:inantrial(i) - 1); NaN; TimeIn_ranking(inantrial(i):end)];
    audStop_ranking = [audStop_ranking(1:2*inantrial(i) - 5); NaN; audStop_ranking(2*inantrial(i) - 4:end)];
    audStop_ranking = [audStop_ranking(1:2*inantrial(i) - 1); NaN; audStop_ranking(2*inantrial(i):end)];
    tacStop_ranking = [tacStop_ranking(1:2*inantrial(i) - 5); NaN; tacStop_ranking(2*inantrial(i) - 4:end)];
    tacStop_ranking = [tacStop_ranking(1:2*inantrial(i) - 1); NaN; tacStop_ranking(2*inantrial(i):end)];
    Choice_P_ranking = [Choice_P_ranking(1:inantrial(i) - 1); NaN; Choice_P_ranking(inantrial(i):end)];
    CAL_ranking = [CAL_ranking(1:inantrial(i) - 1); NaN; CAL_ranking(inantrial(i):end)];
    Head_Change_ranking = [Head_Change_ranking(1:inantrial(i) - 1); NaN; Head_Change_ranking(inantrial(i) + 1:end)];
end

% Reshape
rms_ranked = reshape(rms_ranking,size(rms_total_mat));
TimeIn_ranked = reshape(TimeIn_ranking,size(fracTimeIn_mat));
audStop_ranked = reshape(audStop_ranking,size(audStop_mat));
tacStop_ranked = reshape(tacStop_ranking,size(tacStop_mat));
Choice_P_ranked = reshape(Choice_P_ranking,size(Choice_P_mat));
CAL_ranked = reshape(CAL_ranking,size(CAL_Score_mat));
Head_Change_ranked = reshape(Head_Change_ranking,size(Head_Change_mat));

%% Calculate the Probabilities of each SR level across the distribution
% rms_total_p = zeros(1,length(SRlevels(7:end)));
% fracTimeIn_p = zeros(1,length(SRlevels(7:end)));
% Choice_P_p = zeros(1,length(SRlevels(7:end)));
% CAL_Score_p = zeros(1,length(SRlevels(7:end)));
% Head_Change_p = zeros(1,length(SRlevels(7:end)));
% audStop_p = zeros(1,length(SRlevels(7:end)));
% tacStop_p = zeros(1,length(SRlevels(7:end)));
% 
% for i=1:length(SRlevels(7:end))
%     
%     DenomSize = numel(rms_ranked);
%     for j=1:4
%         rms_total_p(i) = rms_total_p(i)+rms_ranked(j,i)/(DenomSize*4);
%         fracTimeIn_p(i) = fracTimeIn_p(i)+TimeIn_ranked(j,i)/(DenomSize*4);
%         Choice_P_p(i) = Choice_P_p(i)+Choice_P_ranked(j,i)/(DenomSize*4);
%         CAL_Score_p(i) = CAL_Score_p(i)+CAL_ranked(j,i)/(DenomSize*4);
%         Head_Change_p(i) = Head_Change_p(i)+Head_Change_ranked(j,i)/(DenomSize*4);
%     end
%     
%     DenomSize = numel(audStop_ranked);
%     for j=1:8
%         audStop_p(i) = audStop_p(i)+audStop_ranked(j,i)/(DenomSize*8);
%         tacStop_p(i) = tacStop_p(i)+tacStop_ranked(j,i)/(DenomSize*8);
%     end
%     
%     Total_Performance(i) = 1/3*((1-rms_total_p(i)) + (1-fracTimeIn_p(i)) + (1-Head_Change_p(i))) + 1/2*((1-audStop_p(i)) + (1-tacStop_p(i))) + 1/2*(Choice_P_p(i) + CAL_Score_p(i));
% 
% end

rms_total_p = mean(rms_ranked, 'omitnan')/numel(rms_ranked);
fracTimeIn_p = mean(TimeIn_ranked, 'omitnan')/numel(TimeIn_ranked);
Choice_P_p = mean(Choice_P_ranked, 'omitnan')/numel(Choice_P_ranked);
CAL_Score_p = mean(CAL_ranked, 'omitnan')/numel(CAL_ranked);
Head_Change_p = mean(Head_Change_ranked, 'omitnan')/numel(Head_Change_ranked);
audStop_p = mean(audStop_ranked, 'omitnan')/numel(audStop_ranked);
tacStop_p = mean(tacStop_ranked, 'omitnan')/numel(tacStop_ranked);

Total_Performance = 1/3*((1-rms_total_p) + (1-fracTimeIn_p) + (1-Head_Change_p)) + 1/2*((1-audStop_p) + (1-tacStop_p)) + 1/2*(Choice_P_p + CAL_Score_p);

%% Plot Data
size = length(SRlevels(7:end))
XtickMarks = 1:size;
for i = 1:size
    if SRlevels(i+6) == 0
        Xlabel{i} = 0;
    elseif SRlevels(i+6) == 1
        Xlabel{i} = 'MM';
    else
        Xlabel{i} = SRlevels(i+6);
    end
end

figure;
plot(XtickMarks,Total_Performance);
ylabel('Total Performance Score')
xlabel('SR Level') % Change this to condition later
xticks(XtickMarks); xticklabels(Xlabel);

if any(Retest)==0
    
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
else
    addPerData = load(fullfile(LocateFile,['PerformanceData' Subject_str 'Orig.mat']));
    % save all this data
    SRPerfData = {[addPerData.SRPerfData{1} SRlevels(7:end)], [addPerData.SRPerfData{2} rms_total_p], [addPerData.SRPerfData{3} fracTimeIn_p], [addPerData.SRPerfData{4} Choice_P_p], [addPerData.SRPerfData{5} CAL_Score_p], [addPerData.SRPerfData{6} Head_Change_p], [addPerData.SRPerfData{7} audStop_p], [addPerData.SRPerfData{8} tacStop_p], [addPerData.SRPerfData{9} Total_Performance]};
    savefilename = ['PerformanceData' Subject_str 'Retest.mat'];
    saveFilename = fullfile(LocateFile,savefilename);
    save(saveFilename,'SRPerfData');
    
    addRankData = load(fullfile(LocateFile,['RankedData' Subject_str 'Orig.mat']));
    SRRankData = {SRlevels(7:end), rms_ranked, TimeIn_ranked, Choice_P_ranked, CAL_ranked, Head_Change_ranked, audStop_ranked, tacStop_ranked};
    savefilename = ['RankedData' Subject_str 'Retest.mat'];
    saveFilename = fullfile(LocateFile,savefilename);
    save(saveFilename,'SRRankData');    
    
end





