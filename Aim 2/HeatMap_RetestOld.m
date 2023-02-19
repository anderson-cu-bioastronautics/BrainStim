clear all; close all; clc;

%% Get the deets
Subjects = [1 3 5 6];
sizeOrg = 6;
count = 1;
tasks = {'RMS', 'Input', 'Heading', 'Choice', 'CAL', 'Auditory', 'Tactile', 'Total'};
colors = ['r';'b';'m';'g';'c'];

for Subs = Subjects
    
    Subject = Subs;
    currentFolder = pwd;
    
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
    Sortedmeans = mean(RetDataSorted');
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

%% All line plots of all subjects
xlab = {'sham' 'nGVS' 'AWN' 'MMSR'};
xtickmarks = [1:4]; xlimit = [0.5 4.5];
figure()
%colorbar
plotTitle = ['SR Performance'];
sgtitle(plotTitle);

for i = 1:length(tasks)
    subplot(2,length(tasks)/2 ,i);
    for sub = 1:length(Subjects)
        PerformData = cell2mat(AllRetest_mat(sub));
        meansData = cell2mat(means_mat(sub))';
        TaskData(sub,:) = PerformData(i+1,:);
        Taskmean(sub,:) = meansData(i+1,:);
        plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-*' colors(sub)],'LineWidth',1.15)
        hold on
        title(['Task ' tasks{i}]); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);
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
        243 80 80 %
        245 145 145
        240 167 167 %%
        249 195 195
        250 216 216 %
        219 255 218 %
        180 249 195
        159 242 136 %%
        109 229 95
        106 191 82 %
        47 168 34]/255;
    
    yvalues = {'RMS', 'Input', 'Heading', 'Choice', 'CAL', 'Auditory', 'Tactile', 'Total'};
    xvalues = {'nGVS', 'AWN', 'MMSR'};
    
    subplot(1,length(Subjects),sub);
    heatmap(xvalues, yvalues, PercentData)%, 'Colormap',newcolors, 'caxis', ([-60 60]));
    title(['Subject ' num2str(Subjects(sub))]);
    colormap(newcolors)
    caxis([-maxPerc maxPerc])
    
    axs = struct(gca); %ignore warning that this should be avoided
    c = axs.Colorbar;
    c.Ticks = -90:15:90;
    c.TickLabels = {'-90+', '-75', '-60','-45','-30', '-15', '0', '15', '30', '45', '60', '75', '90+'};
    
    
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
        243 80 80 %
        245 145 145
        240 167 167 %%
        249 195 195
        250 216 216 %
        219 255 218 %
        180 249 195
        159 242 136 %%
        109 229 95
        106 191 82 %
        47 168 34]/255;
    
    yvalues = num2cell(num2str(Subjects'));
    xvalues = {'nGVS', 'AWN', 'MMSR'};
    
    subplot(2,length(tasks)/2 ,i);
    heatmap(xvalues, yvalues, TaskData)%, 'Colormap',newcolors, 'caxis', ([-60 60]));
    title(['Task ' tasks{i}]);
    colormap(newcolors)
    caxis([-maxPerc maxPerc])
    
    axs = struct(gca); %ignore warning that this should be avoided
    c = axs.Colorbar;
    c.Ticks = -90:15:90;
    c.TickLabels = {'-90+', '-75', '-60','-45','-30', '-15', '0', '15', '30', '45', '60', '75', '90+'};
    
    PercentData = [];
end