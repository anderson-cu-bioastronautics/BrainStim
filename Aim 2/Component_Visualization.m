%% Data to investigate ranked data and all sorts of stuff

clear all; close all; clc;

%% Get the deets
Subject = input('What is the subject number : ');
Retest = input('Enter \n 0 for Original Test Analysis \n 1 for Retest Analysis \nInput: ');
currentFolder = pwd;

if any(Retest) == 0
    Order = 1:7; xlimited = [0 7];
else
    Order = 1:11; xlimited = [0 11];
end

% Subject file information
if Subject >9
    Subject_str = num2str(Subject);
else
    Subject_str = ['0',num2str(Subject)];
end

load(['PerformanceData' Subject_str '.mat']);
load(['RankedData' Subject_str '.mat']);

Performance_Data = SRPerfData;
Ranked_Data = SRRankData;

%% Parse and Sort the data for SR Analysis
SRlevels = Performance_Data{1,1};
rms_ranked = Ranked_Data{1,2};
TimeIn_ranked = Ranked_Data{1,3};
Choice_P_ranked = Ranked_Data{1,4};
CAL_ranked = Ranked_Data{1,5};
Head_Change_ranked = Ranked_Data{1,6};
audStop_ranked = Ranked_Data{1,7};
tacStop_ranked = Ranked_Data{1,8};

rms_p = Performance_Data{1,2};
TimeIn_p = Performance_Data{1,3};
Choice_P_p = Performance_Data{1,4};
CAL_p = Performance_Data{1,5};
Head_Change_p = Performance_Data{1,6};
audStop_p = Performance_Data{1,7};
tacStop_p = Performance_Data{1,8};
Total_p = Performance_Data{1,9};

%% Plot performance data through time order
h1=figure()
plotTitle = ['Performance in Order - Subject ' Subject_str];
sgtitle(plotTitle);

subplot(2,4,1)
plot(Order,rms_p,'k','LineWidth',1.15)
title('Root Mean Squared'); xlim(xlimited);

subplot(2,4,2)
plot(Order,TimeIn_p,'k','LineWidth',1.15)
title('Throttle Input'); xlim(xlimited);

subplot(2,4,3)
plot(Order,Head_Change_p,'k','LineWidth',1.15)
title('Smoothness'); xlim(xlimited);

subplot(2,4,4)
plot(Order,Choice_P_p,'k','LineWidth',1.15)
title('LP Choice'); xlim(xlimited);

subplot(2,4,5)
plot(Order,CAL_p,'k','LineWidth',1.15)
title('Crash Abort Land'); xlim(xlimited);

subplot(2,4,6)
plot(Order,audStop_p,'k','LineWidth',1.15)
title('Auditory Reaction'); xlim(xlimited);

subplot(2,4,7)
plot(Order,tacStop_p,'k','LineWidth',1.15)
title('Tactile Reaction'); xlim(xlimited);

subplot(2,4,8)
plot(Order,Total_p,'k','LineWidth',1.15)
title('Total Perfomrance Score'); xlim(xlimited);

%% Plot performance data through time order
h1=figure()
plotTitle = ['Performance in Order - Subject ' Subject_str];
sgtitle(plotTitle);

subplot(2,4,1)
plot(Order,rms_p,'k','LineWidth',1.15)
title('Root Mean Squared'); xlim(xlimited);

subplot(2,4,2)
plot(Order,TimeIn_p,'k','LineWidth',1.15)
title('Throttle Input'); xlim(xlimited);

subplot(2,4,3)
plot(Order,Head_Change_p,'k','LineWidth',1.15)
title('Smoothness'); xlim(xlimited);

subplot(2,4,4)
plot(Order,Choice_P_p,'k','LineWidth',1.15)
title('LP Choice'); xlim(xlimited);

subplot(2,4,5)
plot(Order,CAL_p,'k','LineWidth',1.15)
title('Crash Abort Land'); xlim(xlimited);

subplot(2,4,6)
plot(Order,audStop_p,'k','LineWidth',1.15)
title('Auditory Reaction'); xlim(xlimited);

subplot(2,4,7)
plot(Order,tacStop_p,'k','LineWidth',1.15)
title('Tactile Reaction'); xlim(xlimited);

subplot(2,4,8)
plot(Order,Total_p,'k','LineWidth',1.15)
title('Total Perfomrance Score'); xlim(xlimited);

%% Plot Individual Performance Enhancement Data
if (Retest)
Order2 = Order((end-3):end); %retest 
Order = Order(1:(end-4)); %original 

size=length(Order);


RetData = Data(:, (size+1):end);
RetDataNC = DataNC(:, (size+1):end);
RetDataPer = DataPer(:, (size+1):end);

Data = Data(:, 1:size);
DataNC = DataNC(:,1:size);
DataPer = DataPer(:,1:size);

[RetOrder, sorted] = sort(Order2, 'ascend'); %will sort to sham, vsr, mmsr, asr
RetDataSorted = RetData(:, sorted);
RetDataSortedNC = RetDataNC(:, sorted);
RetDataSortedPer = RetDataPer(:, sorted);

end
shamLine=RetDataSorted(1,1).*ones(1,length(vsrLev));
shamLineNC = RetDataSortedNC(1,1).*ones(1,length(vsrLev));
shamLinePer = RetDataSortedPer(1,1).*ones(1,length(vsrLev));
