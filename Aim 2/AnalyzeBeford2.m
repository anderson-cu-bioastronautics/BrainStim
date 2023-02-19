clear all; close all; clc;

%% Read in .mat files
%each row is a subject, each column is nGVS, AWN, MMSR, each panel is for a test (DSST, PVT, etc)

Subjects = [10:17 26:27 29:33];
numSubs = length(Subjects);

currentFolder = pwd;
folder = [pwd, '\', 'Subject_Data'];
saveFolder = [pwd, '\', 'plots'];
mkdir(saveFolder);

for sub = 1:length(Subjects)
    files = fullfile(folder, ['Subject_' num2str(Subjects(sub))]);
    bedData = xlsread(fullfile(files,['Subject_' num2str(Subjects(sub)) '_Bedford_Data.xls']));
    
    load(fullfile(files,['PerformanceData' num2str(Subjects(sub)) 'Retest.mat']));
    levelSub = SRPerfData{1,1}(7:10);
    [levSort,SortInd] = sort(levelSub); % Need to sort data in an appropriate way
    SortWant = [SortInd(1:2) SortInd(4) SortInd(3)];
    levels = levelSub(SortWant);
    BedWant = bedData(7:10,2);
    BedDWant = BedWant(SortWant);
    
    SubLevels(sub,:) = levels;
    BedMat(sub,:) = BedDWant;
    Cond_numMat(sub,:) = [1:4];
end

% Organize data for RMANOVA later
Sub_vect = 1:numSubs; Sub_vect = Sub_vect';
Sub_mat = [Sub_vect Sub_vect Sub_vect Sub_vect];
Condition = reshape(Cond_numMat,[numSubs*4,1]);
Subject = reshape(Sub_mat,[numSubs*4,1]);
Bedford = reshape(BedMat,[numSubs*4,1]);

%% Plot Bedford Data
tickL = 1:4;
markers = ['o' '+' '*' '.' 'x' 'p' 'h' '^' '>' '<' 'v' 's' 'd' 'o' '+'];
alpha = 0.4;
xvalues = {'Sham', 'nGVS', 'AWN', 'MMSR'};

figure()
plotTitle = ['Bedford Performance across Treatments'];
sgtitle(plotTitle);

% for sub = 1:numSubs
%     plot(tickL, BedMat(sub,:),'*','color',[0,0,0]+alpha, 'Marker', markers(sub))
%     hold on
% end

yMean = mean(BedMat(:,1:4));
ySEM = std(BedMat(:,1:4))/sqrt(numSubs);
CI95 = tinv([0.025 0.975], numSubs-1);
yCI95 = bsxfun(@times, ySEM, CI95(:));
plot(tickL, yMean ,'-k','LineWidth',1.35)
errorbar(tickL, yMean, yCI95(1,:), yCI95(2,:) ,'-k','LineWidth',1.35)
hold off
ylabel('Bedford Score')
xticks(tickL); xticklabels(xvalues);
xlim([0.75 4.25]); ylim([1 5]);
set(gca, 'YTick', 1:5)
set(gca,'FontSize',13)

figName = ['GrayPlot_Bedford.jpg'];
saveas(gcf, fullfile(saveFolder, figName));

%% Run statistical analysis on Bedford data

[p,tbl,stats,terms] = anovan(Bedford,{Condition,Subject},'random',2,'varnames',{'Condition','Subject'});
assumptions_checkingBed(p,tbl,stats,terms,Bedford,Condition,Subject)

%% Run Friedman, cause why not
[pF, tblF] = friedman(BedMat,1)