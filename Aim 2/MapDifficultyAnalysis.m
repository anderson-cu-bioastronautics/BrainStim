clear all; close all; clc;

%% Get the deets
Subjects = [10:17 26:33];
sizeOrg = 6;
count = 1;
tasks = {'RMS', 'Stick Input', 'Smoothness', 'LP Choice', 'Crash Abort Land', 'Auditory Reaction', 'Tactile Reaction', 'Total Score'};;
% colors = ['r';'b';'m';'g';'c'];
colors = linspace(0, 0.75, length(Subjects))'*ones(1, 3);
symbs = 'oxds+*^vph<oxds+*^vph<';
MapDiffInd = [1 3 2 2; 1 2 2 3; 2 3 1 2; 2 2 3 1]'; % Map difficulty indices
[MapSortVal,MapSortInd] = sort(MapDiffInd,1);
PerMapDiffInd = [1 3 2 2 1 3 2 2; 1 2 2 3 1 2 2 3; 2 3 1 2 2 3 1 2; 2 2 3 1 2 2 3 1]'; % Map difficulty indices
[PerMapSortVal,PerMapSortInd] = sort(PerMapDiffInd,1);


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
    
    load(['Subject_Data/Subject_' Subject_str '/IndividualPerformanceData' Subject_str 'Retest.mat']);
    
    Performance_DataLev = SRPerfData;    
    Performance_Data = SRPerfIndData;
    %Ranked_Data = SRRankData;
    
    SRlev = Performance_DataLev{1,1};
    rms_p1 = 1-Performance_Data{1,1};
    TimeIn_p1 = 1-Performance_Data{1,2};
    Choice_P_p1 = Performance_Data{1,6};
    CAL_p1 = Performance_Data{1,7};
    Head_Change_p1 = 1-Performance_Data{1,3};
    audStop_p1 = 1-Performance_Data{1,4};
    tacStop_p1 = 1-Performance_Data{1,5};

    for i = 1:4
    rms_p2(:,i) = rms_p1(MapSortInd(:,i),i);
    TimeIn_p2(:,i) = TimeIn_p1(MapSortInd(:,i),i);
    Choice_P_p2(:,i) = Choice_P_p1(MapSortInd(:,i),i);
    CAL_p2(:,i) = CAL_p1(MapSortInd(:,i),i);
    Head_Change_p2(:,i) = Head_Change_p1(MapSortInd(:,i),i);
    audStop_p2(:,i) = audStop_p1(PerMapSortInd(:,i),i);
    tacStop_p2(:,i) = tacStop_p1(PerMapSortInd(:,i),i);
    end
    

    SRlevels{count} = SRlev((end-3):end);
    rms_p{count} = rms_p2;
    TimeIn_p{count} = TimeIn_p2;
    Choice_P_p{count} = Choice_P_p2;
    CAL_p{count} = CAL_p2;
    Head_Change_p{count} = Head_Change_p2;
    audStop_p{count} = audStop_p2;
    tacStop_p{count} = tacStop_p2;
    
    
    %Total_p = Performance_Data{1,9};
    PerfData{count} = {SRlevels{count}; rms_p{count}; TimeIn_p{count};  Head_Change_p{count}; Choice_P_p{count}; CAL_p{count}; audStop_p{count}; tacStop_p{count}};
    
    Order = SRlevels{count};
    Order2 = Order((end-3):end); %retest
    
    RetData{count} = PerfData{count};
      
    [RetOrder, sorted] = sort(Order2, 'ascend'); %will sort to sham, vsr, mmsr, asr
    sorted(:,3:4) = flip(sorted(:,3:4),2);
    RetDataL = RetData{1,count};
    
    for i = 1:8
    RetDataSorted{i} = RetDataL{i}(:, sorted);
    Sortedmeans{i} = mean(RetDataSorted{i}','omitnan');
    end
    AllRetest_mat{count} = RetDataSorted;
   % AllPercent_mat{count} = PercentData;
    means_mat{count} = Sortedmeans;
    PercentData = []; RetDataL = []; PerfData = []; Performance_Data =[]; RetDataSorted = [];
    SRlev = []; rms_p2 = []; TimeIn_p2 = []; Choice_P_p2 = []; CAL_p2 = []; Head_Change_p2 = []; audStop_p2 = []; tacStop_p2 = [];
    count = count+1;
end

%% Easy line plots of all subjects (symbol coding)
xlab = {'sham' 'nGVS' 'AWN' 'MMSR'};
xtickmarks = [1:4]; xlimit = [0.5 4.5];
figure()
set(gca,'FontSize',14)
%colorbar
plotTitle = ['SR Performance on Easy Maps'];
sgtitle(plotTitle);

for i = 1:7
    subplot(2,4 ,i);
    for sub = 1:8
        PerformData = (AllRetest_mat{sub}{i+1}); % subject then task
        meansData = (means_mat{sub}{i+1})';
        if i < 6
        TaskData(sub,:) = PerformData(1,:);
        Taskmean(sub,:) = meansData(1,:);
        else
        TaskData(sub,:) = mean(PerformData(1:2,:),'omitnan');
        Taskmean(sub,:) = mean(meansData(1:2,:),'omitnan');            
        end
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-*' colors(sub)],'LineWidth',1.15)
        plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),[ symbs(sub) 'k'],'LineWidth',1.15, 'MarkerFaceColor', 'k')
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),'-*','LineWidth',1.15, 'Color', colors(sub, :), 'MarkerFaceColor', colors(sub, :))
        hold on
        title(tasks{i}); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);
    end
        plot(xtickmarks,mean(TaskData(1:8,:),'omitnan')-mean(Taskmean(1:8,:),'omitnan'),[ '-k'],'LineWidth',1.25)
    for sub = 9:16
        PerformData = (AllRetest_mat{sub}{i+1}); % subject then task
        meansData = (means_mat{sub}{i+1})';
        if i < 6
        TaskData(sub,:) = PerformData(1,:);
        Taskmean(sub,:) = meansData(1,:);
        else
        TaskData(sub,:) = mean(PerformData(1:2,:),'omitnan');
        Taskmean(sub,:) = mean(meansData(1:2,:),'omitnan');            
        end
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-*' colors(sub)],'LineWidth',1.15)
        plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),[ symbs(sub) 'r'],'LineWidth',1.15, 'MarkerFaceColor', 'r')
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),'-*','LineWidth',1.15, 'Color', colors(sub, :), 'MarkerFaceColor', colors(sub, :))
        hold on
        title(tasks{i}); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);
    end
        plot(xtickmarks,mean(TaskData(9:16,:),'omitnan')-mean(Taskmean(9:16,:),'omitnan'),[ '-r'],'LineWidth',1.25)    
    yline(0,'LineWidth',.4)
    hold off
    set(gca,'FontSize',13)
    TaskData(any(isnan(TaskData), 2), :) = [];
    Taskmean(any(isnan(Taskmean), 2), :) = [];
    EasyData{i} = TaskData; % by task
    EasyDataMean{i} = Taskmean; % by task 
    if i <6 % not auditory and tactile measurements
    [Easyp{i}, tbl, Easystats{i}] = friedman(TaskData,1,'off');
    else
    [Easyp{i}, tbl, Easystats{i}] = friedman(TaskData,1,'off');
    end
    TaskData = [];
        Taskmean = [];    
end

%% Hard line plots of all subjects (symbol coding)
xlab = {'sham' 'nGVS' 'AWN' 'MMSR'};
xtickmarks = [1:4]; xlimit = [0.5 4.5];
figure()
set(gca,'FontSize',14)
%colorbar
plotTitle = ['SR Performance on Hard Maps'];
sgtitle(plotTitle);

for i = 1:7
    subplot(2,4 ,i);
    for sub = 1:8
        PerformData = (AllRetest_mat{sub}{i+1}); % subject then task
        meansData = (means_mat{sub}{i+1})';
        if i < 6
        TaskData(sub,:) = PerformData(4,:);
        Taskmean(sub,:) = meansData(4,:);
        else
        TaskData(sub,:) = mean(PerformData(7:8,:),'omitnan');
        Taskmean(sub,:) = mean(meansData(7:8,:),'omitnan');            
        end
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-*' colors(sub)],'LineWidth',1.15)
        plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),[ symbs(sub) 'k'],'LineWidth',1.15, 'MarkerFaceColor', 'k')
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),'-*','LineWidth',1.15, 'Color', colors(sub, :), 'MarkerFaceColor', colors(sub, :))
        hold on
        title(tasks{i}); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);
    end
        plot(xtickmarks,mean(TaskData(1:8,:),'omitnan')-mean(Taskmean(1:8,:),'omitnan'),[ '-k'],'LineWidth',1.25)
    for sub = 9:16
        PerformData = (AllRetest_mat{sub}{i+1}); % subject then task
        meansData = (means_mat{sub}{i+1})';
        if i < 6
        TaskData(sub,:) = PerformData(4,:);
        Taskmean(sub,:) = meansData(4,:);
        else
        TaskData(sub,:) = mean(PerformData(7:8,:),'omitnan');
        Taskmean(sub,:) = mean(meansData(7:8,:),'omitnan');            
        end
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-*' colors(sub)],'LineWidth',1.15)
        plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),[ symbs(sub) 'r'],'LineWidth',1.15, 'MarkerFaceColor', 'r')
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),'-*','LineWidth',1.15, 'Color', colors(sub, :), 'MarkerFaceColor', colors(sub, :))
        hold on
        title(tasks{i}); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);
    end
        plot(xtickmarks,mean(TaskData(9:16,:),'omitnan')-mean(Taskmean(9:16,:),'omitnan'),[ '-r'],'LineWidth',1.25)    
    yline(0,'LineWidth',.4)
    hold off
    set(gca,'FontSize',13)
    TaskData(any(isnan(TaskData), 2), :) = [];
    Taskmean(any(isnan(Taskmean), 2), :) = [];
    HardData{i} = TaskData; % by task
    HardDataMean{i} = Taskmean; % by task 
    if i <6 % not auditory and tactile measurements
    [Hardp{i}, tbl, Hardstats{i}] = friedman(TaskData,1,'off');
    else
    [Hardp{i}, tbl, Hardstats{i}] = friedman(TaskData,1,'off');
    end
        TaskData = [];
        Taskmean = [];    
end

%% Medium line plots of all subjects (symbol coding)
xlab = {'sham' 'nGVS' 'AWN' 'MMSR'};
xtickmarks = [1:4]; xlimit = [0.5 4.5];
figure()
set(gca,'FontSize',14)
%colorbar
plotTitle = ['SR Performance on Medium Maps'];
sgtitle(plotTitle);

for i = 1:7
    subplot(2,4 ,i);
    for sub = 1:8
        PerformData = (AllRetest_mat{sub}{i+1}); % subject then task
        meansData = (means_mat{sub}{i+1})';
        if i < 6
        TaskData(sub,:) = mean(PerformData(2:3,:),'omitnan');
        Taskmean(sub,:) = mean(meansData(2:3,:),'omitnan');
        else
        TaskData(sub,:) = mean(PerformData(3:6,:),'omitnan');
        Taskmean(sub,:) = mean(meansData(3:6,:),'omitnan');            
        end
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-*' colors(sub)],'LineWidth',1.15)
        plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),[ symbs(sub) 'k'],'LineWidth',1.15, 'MarkerFaceColor', 'k')
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),'-*','LineWidth',1.15, 'Color', colors(sub, :), 'MarkerFaceColor', colors(sub, :))
        hold on
        title(tasks{i}); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);
    end
        plot(xtickmarks,mean(TaskData(1:8,:),'omitnan')-mean(Taskmean(1:8,:),'omitnan'),[ '-k'],'LineWidth',1.25)
    for sub = 9:16
        PerformData = (AllRetest_mat{sub}{i+1}); % subject then task
        meansData = (means_mat{sub}{i+1})';
        if i < 6
        TaskData(sub,:) = mean(PerformData(2:3,:),'omitnan');
        Taskmean(sub,:) = mean(meansData(2:3,:),'omitnan');
        else
        TaskData(sub,:) = mean(PerformData(3:6,:),'omitnan');
        Taskmean(sub,:) = mean(meansData(3:6,:),'omitnan');            
        end
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),['-*' colors(sub)],'LineWidth',1.15)
        plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),[ symbs(sub) 'r'],'LineWidth',1.15, 'MarkerFaceColor', 'r')
%         plot(xtickmarks,TaskData(sub,:)-Taskmean(sub),'-*','LineWidth',1.15, 'Color', colors(sub, :), 'MarkerFaceColor', colors(sub, :))
        hold on
        title(tasks{i}); xticks(xtickmarks); xticklabels(xlab); xlim(xlimit);
    end
        plot(xtickmarks,mean(TaskData(9:16,:),'omitnan')-mean(Taskmean(9:16,:),'omitnan'),[ '-r'],'LineWidth',1.25)    
    yline(0,'LineWidth',.4)
    hold off
    set(gca,'FontSize',13)
    TaskData(any(isnan(TaskData), 2), :) = [];
    Taskmean(any(isnan(Taskmean), 2), :) = [];
    MediumData{i} = TaskData; % by task
    MediumDataMean{i} = Taskmean; % by task 
    if i <6 % not auditory and tactile measurements
    [Mediump{i}, tbl, Mediumstats{i}] = friedman(TaskData,1,'off');
    else
    [Mediump{i}, tbl, Mediumstats{i}] = friedman(TaskData,1,'off');
    end 
        TaskData = [];
        Taskmean = [];    
end