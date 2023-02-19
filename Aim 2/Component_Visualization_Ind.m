%% Data to investigate ranked data and all sorts of stuff

clear all; close all;

%% Get the deets
Subject = input('What is the subject number : ');
Retest = input('Enter \n 0 for Original Test Analysis \n 1 for Retest Analysis \nInput: ');
currentFolder = pwd;

sizeOrg = 6;

if any(Retest) == 0
    Order = 1:6; xlimited = [0.5 6.5];
else
    Order = 1:10; xlimited = [0.5 10.5];
end

% Subject file information
if Subject >9
    Subject_str = num2str(Subject);
else
    Subject_str = ['0',num2str(Subject)];
end

load(['PerformanceData' Subject_str '.mat']);
%load(['RankedData' Subject_str '.mat']);

Performance_Data = SRPerfData;
%Ranked_Data = SRRankData;

%% Parse and Sort the data for SR Analysis
SRlevels = Performance_Data{1,1};
% rms_ranked = Ranked_Data{1,2};
% TimeIn_ranked = Ranked_Data{1,3};
% Choice_P_ranked = Ranked_Data{1,4};
% CAL_ranked = Ranked_Data{1,5};
% Head_Change_ranked = Ranked_Data{1,6};
% audStop_ranked = Ranked_Data{1,7};
% tacStop_ranked = Ranked_Data{1,8};

rms_p = 1-Performance_Data{1,2};
TimeIn_p = 1-Performance_Data{1,3};
Choice_P_p = Performance_Data{1,4};
CAL_p = Performance_Data{1,5};
Head_Change_p = 1-Performance_Data{1,6};
audStop_p = 1-Performance_Data{1,7};
tacStop_p = 1-Performance_Data{1,8};
Total_p = Performance_Data{1,9};
PerfData = [SRlevels; rms_p; TimeIn_p;  Head_Change_p; Choice_P_p; CAL_p; audStop_p; tacStop_p; Total_p]

%% Plot performance data through time order

for i = 1:length(Order)
    if Order(i) == 0
        xlab{i} = 0;
    elseif PerfData(1,i) == 1
        xlab{i} = 'MM';
    else
        xlab{i} = PerfData(1,i);
    end
end

xtickmarks = [Order];

h1=figure();
plotTitle = ['Performance in Order - Subject ' Subject_str];
sgtitle(plotTitle);

subplot(2,4,1)
plot(Order,rms_p,'k','LineWidth',1.15)
title('Root Mean Squared'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimited);

subplot(2,4,2)
plot(Order,TimeIn_p,'k','LineWidth',1.15)
title('Stick Input'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimited);

subplot(2,4,3)
plot(Order,Head_Change_p,'k','LineWidth',1.15)
title('Smoothness'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimited);

subplot(2,4,4)
plot(Order,Choice_P_p,'k','LineWidth',1.15)
title('LP Choice'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimited);

subplot(2,4,5)
plot(Order,CAL_p,'k','LineWidth',1.15)
title('Crash Abort Land');xticks(xtickmarks); xticklabels(xlab); xlim(xlimited);

subplot(2,4,6)
plot(Order,audStop_p,'k','LineWidth',1.15)
title('Auditory Reaction'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimited);

subplot(2,4,7)
plot(Order,tacStop_p,'k','LineWidth',1.15)
title('Tactile Reaction');xticks(xtickmarks); xticklabels(xlab);  xlim(xlimited);

subplot(2,4,8)
plot(Order,Total_p,'k','LineWidth',1.15)
title('Total Performance Score'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimited);


%% Plot Individual Performance Enhancement Data for Retests Only
if any(Retest) == 1

Order = SRlevels;      
Order2 = Order((end-3):end); %retest 

RetData = PerfData(:, (sizeOrg+1):end);

OrigData = PerfData(:, 1:sizeOrg);

[RetOrder, sorted] = sort(Order2, 'ascend'); %will sort to sham, vsr, mmsr, asr
RetDataSorted = RetData(:, sorted);
RetDataSorted(:,3:4) = flip(RetDataSorted(:,3:4),2);
Sortedmeans = mean(RetDataSorted');


xlab = {'sham' RetOrder(2) RetOrder(4) 'MMSR'};
xtickmarks = [1:4]; xlimit = [0.5 4.5];

h1=figure();
plotTitle = ['SR Performance - Subject ' Subject_str];
sgtitle(plotTitle);

subplot(2,4,1)
plot(xtickmarks,RetDataSorted(2,:)-Sortedmeans(2),'-*k','LineWidth',1.15)
title('Root Mean Squared'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit); 

subplot(2,4,2)
plot(xtickmarks,RetDataSorted(3,:)-Sortedmeans(3),'-*k','LineWidth',1.15)
title('Stick Input'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);

subplot(2,4,3)
plot(xtickmarks,RetDataSorted(4,:)-Sortedmeans(4),'-*k','LineWidth',1.15)
title('Smoothness'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);

subplot(2,4,4)
plot(xtickmarks,RetDataSorted(5,:)-Sortedmeans(5),'-*k','LineWidth',1.15)
title('LP Choice'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);

subplot(2,4,5)
plot(xtickmarks,RetDataSorted(6,:)-Sortedmeans(6),'-*k','LineWidth',1.15)
title('Crash Abort Land'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);

subplot(2,4,6)
plot(xtickmarks,RetDataSorted(7,:)-Sortedmeans(7),'-*k','LineWidth',1.15)
title('Auditory Reaction'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);

subplot(2,4,7)
plot(xtickmarks,RetDataSorted(8,:)-Sortedmeans(8),'-*k','LineWidth',1.15)
title('Tactile Reaction'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);

subplot(2,4,8)
plot(xtickmarks,RetDataSorted(9,:)-Sortedmeans(9),'-*k','LineWidth',1.15)
title('Total Performance Score'); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);

end