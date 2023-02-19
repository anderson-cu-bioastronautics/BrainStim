clear all;
% close all; clc;

%% Get the deets
% Subjects = [26:33];
% Subjects = [10:17];
Subjects = [10:17 26:33];
sizeOrg = 6;
count = 1;
tasks = {'RMS', 'Stick Input', 'Smoothness', 'LP Choice', 'Crash Abort Land', 'Auditory Reaction', 'Tactile Reaction', 'Total Score'};
% colors = ['r';'b';'m';'g';'c'];
colors = linspace(0, 0.75, length(Subjects))'*ones(1, 3);
symbs = 'oxds+*^vph<';

for Subs = Subjects
    
    Subject = Subs;
    currentFolder = pwd;
    
    % Subject file information
    if Subject >9
        Subject_str = num2str(Subject);
    else
        Subject_str = ['0',num2str(Subject)];
    end
    
    load(['Subject_Data/Subject_' Subject_str '/PerformanceData' Subject_str 'Retest.mat']);
    %load(['RankedData' Subject_str '.mat']);
    
    Performance_Data = SRPerfData;
    %Ranked_Data = SRRankData;
    
    SRlevels = Performance_Data{1,1};
    rms_p = 1-Performance_Data{1,2};
    TimeIn_p = 1-Performance_Data{1,3};
    Choice_P_p = Performance_Data{1,4};
    CAL_p = Performance_Data{1,5};
    Head_Change_p = 1-Performance_Data{1,6};
    audStop_p = 1-Performance_Data{1,7};
    tacStop_p = 1-Performance_Data{1,8};
    Total_p = Performance_Data{1,9};
    PerfData = [SRlevels; rms_p; TimeIn_p;  Head_Change_p; Choice_P_p; CAL_p; audStop_p; tacStop_p; Total_p];
    
    Order = SRlevels;
    Order2 = Order((end-3):end); %retest
    
    RetData = PerfData(:, (sizeOrg+1):end);
    
    OrigData = PerfData(:, 1:sizeOrg);
    
    [RetOrder, sorted] = sort(Order2, 'ascend'); %will sort to sham, vsr, mmsr, asr
    RetDataSorted = RetData(:, sorted);
    Sortedmeans = mean(RetDataSorted', 1, 'omitnan');
    RetDataSorted(:,3:4) = flip(RetDataSorted(:,3:4),2);
    
    for j=2:9
        for i = 2:4
            PercentData(j-1,i-1) = (RetDataSorted(j,i)-RetDataSorted(j,1))/RetDataSorted(j,1)*100;
        end
    end
    
    AllRetest_mat{count} = RetDataSorted;
    AllPercent_mat{count} = PercentData;
    means_mat{count} = Sortedmeans;
    PercentData = []; RetDataSorted = []; PerfData = []; Performance_Data =[];
    count = count+1;
end

%% All line plots of all subjects (symbol coding)
xlab = {'sham' 'nGVS' 'AWN' 'MMSR'};
xtickmarks = [1:4]; xlimit = [0.5 4.5];
figure()
%colorbar
plotTitle = ['SR Performance'];
sgtitle(plotTitle);

for i = 1:length(tasks)
    subplot(2,length(tasks)/2 ,i);
    subNAA = 1;
    subAA = 1;
    for sub = 1:length(Subjects)
        PerformData = cell2mat(AllRetest_mat(sub));
        meansData = cell2mat(means_mat(sub))';
        TaskData(sub,:) = PerformData(i+1,:);
        Taskmean(sub,:) = meansData(i+1,:);
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-*' colors(sub)],'LineWidth',1.15)
        if Subjects(sub) < 20
            plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-' symbs(subNAA) 'k'],'LineWidth',1.15, 'MarkerFaceColor', 'k', 'MarkerSize', 10);
            subNAA = subNAA + 1;
        else
            plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-' symbs(subAA) 'r'],'LineWidth',1.15, 'MarkerFaceColor', 'r', 'MarkerSize', 10);
            subAA = subAA + 1;
        end
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),'-*','LineWidth',1.15, 'Color', colors(sub, :), 'MarkerFaceColor', colors(sub, :))
        hold on
        title(tasks{i}); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);
    end
    yline(0,'LineWidth',.4)
    hold off
        TaskData = [];
        Taskmean = [];    
end

%% Heat maps per subject
figure()
%colorbar
plotTitle = ['% Difference Optimal Conditions and Sham by Subject '];
sgtitle(plotTitle);

for sub = 1:length(Subjects)
    
    PercentData = cell2mat(AllPercent_mat(sub));
    maxPerc = 90;%max([abs(SpdPerc) abs(AccPerc) abs(FeedPerc)],[], 'all');
    
newcolors = [185 26 26
220 40 40 %%%
243 80 80 %
245 145 145
242 155 155 %%%
240 167 167 %%
249 195 195
250 216 216 %
255 255 255
219 255 218 %
180 249 195
159 242 136 %%
130 235 120 %%%
109 229 95
106 191 82 %
55 180 50 %%%
47 168 34]/255;
    
    yvalues = {'RMS', 'Stick Input', 'Smoothness', 'LP Choice', 'Crash Abort Land', 'Auditory Reaction', 'Tactile Reaction', 'Total Score'};
    xvalues = {'nGVS', 'AWN', 'MMSR'};
    
    subplot(1,length(Subjects),sub);
    heatmap(xvalues, yvalues, PercentData, 'CellLabelFormat', '%0.1f')%, 'Colormap',newcolors, 'caxis', ([-60 60]));
    if sub > 1
        set(gca, 'YDisplayLabels', {'', '', '', '', '', '', '', ''});
    end
    title(['Subject ' num2str(Subjects(sub))]);
    colormap(newcolors)
    caxis([-maxPerc maxPerc])
    
    axs = struct(gca); %ignore warning that this should be avoided
    c = axs.Colorbar;
    c.Ticks = -85:10:85;
    c.TickLabels = {'-85+','-75', '-65', '-55','-45','-35', '-25', '-15', '-5', '5', '15', '25', '35', '45', '55', '65', '75','85+'};
    
    
    PercentData = [];
end

%% Heat maps per Task
figure()
%colorbar
plotTitle = ['% Difference Optimal Conditions and Sham by Task '];
sgtitle(plotTitle);

for i = 1:length(tasks)
    
    for sub = 1:length(Subjects)
        PercentData = cell2mat(AllPercent_mat(sub));
        TaskData(sub,:) = PercentData(i,:);
    end
    
    maxPerc = 90;%max([abs(SpdPerc) abs(AccPerc) abs(FeedPerc)],[], 'all');
    
newcolors = [185 26 26
220 40 40 %%%
243 80 80 %
245 145 145
242 155 155 %%%
240 167 167 %%
249 195 195
250 216 216 %
255 255 255
219 255 218 %
180 249 195
159 242 136 %%
130 235 120 %%%
109 229 95
106 191 82 %
55 180 50 %%%
47 168 34]/255;
%FOR -85:10:85
    
    yvalues = num2cell(Subjects');
    xvalues = {'nGVS', 'AWN', 'MMSR'};
    
    subplot(2,length(tasks)/2 ,i);
    heatmap(xvalues, yvalues, TaskData, 'CellLabelFormat', '%0.1f')%, 'Colormap',newcolors, 'caxis', ([-60 60]));
    title(tasks{i});
    colormap(newcolors)
    caxis([-maxPerc maxPerc])
    
    axs = struct(gca); %ignore warning that this should be avoided
    c = axs.Colorbar;
    c.Ticks = -85:10:85;
    c.TickLabels = {'-85+','-75', '-65', '-55','-45','-35', '-25', '-15', '-5', '5', '15', '25', '35', '45', '55', '65', '75','85+'};
    
    PercentData = [];
end