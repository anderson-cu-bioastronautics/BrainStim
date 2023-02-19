%% ANOVA for Metrics Across All Maps
% Y. Shen
% 18 August 2021

%% Load data
clear all; clc; %delete(findall(0));
SubjectsNAA = 10:17;
SubjectsAA = 26:33;
Subjects = [SubjectsNAA SubjectsAA];
sizeOrg = 6;
metrics = {'RMS', 'Stick Input', 'Smoothness', 'LP Choice', 'Crash Abort Land', 'Auditory Reaction', 'Tactile Reaction', 'Total Score'};
N_metrics = 8;
conds = {'sham' 'nGVS' 'AWN' 'MMSR'};
N_conds = 4;
countNAA = 1;
countAA = 1;

AllRetest_mat_NAA = zeros(length(SubjectsNAA), 9, 4);
means_mat_NAA = zeros(length(SubjectsNAA), 9, 1);
AllRetest_mat_AA = zeros(length(SubjectsAA), 9, 4);
means_mat_AA = zeros(length(SubjectsAA), 9, 1);

for Subject = [SubjectsNAA SubjectsAA]
    
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
    
    %     for j=2:9
    %         for i = 2:4
    %             PercentData(j-1,i-1) = (RetDataSorted(j,i)-RetDataSorted(j,1))/RetDataSorted(j,1)*100;
    %         end
    %     end
    if Subject < 20
        AllRetest_mat_NAA(countNAA, :, :) = RetDataSorted;
        %     AllPercent_mat{count} = PercentData;
        means_mat_NAA(countNAA, :, :) = Sortedmeans;
        countNAA = countNAA + 1;
    else
        AllRetest_mat_AA(countAA, :, :) = RetDataSorted;
        %     AllPercent_mat{count} = PercentData;
        means_mat_AA(countAA, :, :) = Sortedmeans;
        countAA = countAA + 1;
    end
    PercentData = []; RetDataSorted = []; PerfData = []; Performance_Data =[];
end

% Drop SR condition row
AllRetest_mat_NAA(:, 1, :) = [];
means_mat_NAA(:, 1, :) = [];
AllRetest_mat_AA(:, 1, :) = [];
means_mat_AA(:, 1, :) = [];

% Deduct means
AllRetest_Centered_NAA = AllRetest_mat_NAA - means_mat_NAA;
AllRetest_Centered_AA = AllRetest_mat_AA - means_mat_AA;

%% Check ANOVA Assumptions

% Test for normality

fprintf('\nNon-Astronaut-Aged:\n');
for i = 1:N_metrics
    xn = reshape((AllRetest_Centered_NAA(:, i, :) - mean(AllRetest_Centered_NAA(:, i, :), 'all'))/std(AllRetest_Centered_NAA(:, i, :), 0, 'all'), length(SubjectsNAA)*N_conds, 1);
    [h, p] = kstest(xn);
    fprintf('%s - H = %g, p = %g\n', metrics{i}, h, p);
end

fprintf('\nAstronaut-Aged:\n');
for i = 1:N_metrics
    xn = reshape((AllRetest_Centered_AA(:, i, :) - mean(AllRetest_Centered_AA(:, i, :), 'all'))/std(AllRetest_Centered_AA(:, i, :), 0, 'all'), length(SubjectsAA)*N_conds, 1);
    [h, p] = kstest(xn);
    fprintf('%s - H = %g, p = %g\n', metrics{i}, h, p);
end

% Check homoscedasticity

varNAA = zeros(N_metrics, N_conds);
varAA = zeros(N_metrics, N_conds);

for i = 1:N_metrics
    for j = 1:N_conds
        varNAA(i, j) = var(AllRetest_Centered_NAA(:, i, j));
        varAA(i, j) = var(AllRetest_Centered_AA(:, i, j));
    end
end

figure;
hold on;
for i = 1:N_conds
    plot(1:N_metrics, varNAA(:, i), 'kd', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    plot(1:N_metrics, varAA(:, i), 'kd', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
end
hold off;
xticklabels(metrics);
xlabel('Metric');
ylabel('Variance');

%% ANOVAs

AllRetest_Centered = [AllRetest_Centered_NAA; AllRetest_Centered_AA];

agegroup = cell(length(Subjects), 1);
for i = 1:length(Subjects)
    if Subjects(i) < 20
        agegroup{i} = 'NAA';
    else
        agegroup{i} = 'AA';
    end
end

fprintf('\nRMANOVAs\n');
fprintf('==========\n');
rm = cell(N_metrics, 1);
for i = 1:N_metrics
    t = table(agegroup, AllRetest_Centered(:, i, 1), AllRetest_Centered(:, i, 2), AllRetest_Centered(:, i, 3), AllRetest_Centered(:, i, 4),...
        'VariableNames', {'agegroup', 'sham', 'nGVS', 'AWN', 'MMSR'});
    condtable = table([1 2 3 4]', 'VariableNames', {'SR Condition'});
    
    fprintf('%s:\n', metrics{i});
%     rm{i} = fitrm(t, 'sham-MMSR~agegroup', 'WithinDesign', condtable);
%     ranova(rm{i})
    
    Performance_mat = t{:,2:5};
    for j=1:length(Performance_mat(:,1))
        Cond_mat(j,:) = [1 2 3 4];
        Sub_mat(j,:) = ones(1,4)*Subjects(j);
        if j <= length(SubjectsNAA)
            Age_mat(j,:) = ones(1,4)*1;
        else
            Age_mat(j,:) = ones(1,4)*2;
        end
    end
    Performance = reshape(Performance_mat,length(Performance_mat(:,1))*4,1);
    Condition = reshape(Cond_mat,length(Performance_mat(:,1))*4,1);
    Age = reshape(Age_mat,length(Performance_mat(:,1))*4,1);
    Subject = reshape(Sub_mat,length(Performance_mat(:,1))*4,1);
    
    tblDat = table(Performance,Condition,Subject);
    if i == 1
    writetable(tblDat,'AvgDataRMS.csv')
    elseif i == 2
    writetable(tblDat,'AvgDataTimeIn.csv')
    elseif i == 3
    writetable(tblDat,'AvgDataSmooth.csv')
    elseif i == 4
    writetable(tblDat,'AvgDataLPC.csv')
    elseif i == 5
    writetable(tblDat,'AvgDataCAL.csv')
    elseif i == 6
    writetable(tblDat,'AvgDataAud.csv')
    elseif i == 7
    writetable(tblDat,'AvgDataTact.csv')
    elseif i == 8
    writetable(tblDat,'AvgDataComp.csv')
    end
    
    %[pA,tblA,statsA,termsA] = anovan(Performance,{Condition,Age,Subject},'model',[1,0,0; 0,1,0; 0,0,1],'nested',[0,0,0; 0,0,0; 0,1,0],'random',3,'varnames',{'Condition','Age','Subject'});
    if i==8
    [pA,tblA,statsA,termsA] = anovan(Performance,{Condition,Age,Subject},'model',[1,0,0; 0,1,0; 0,0,1],'nested',[0,0,0; 0,0,0; 0,1,0],'random',3,'varnames',{'Condition','Age','Subject'});
    end
    Condition = categorical(Condition);
    Subject = categorical(Subject);    
    [pM{i},tblM{i},statsM{i},termsM{i}] = anovan(Performance,{Condition,Subject},'model',[1 0; 0 1 ; 1 1 ],'random',2,'varnames',{'Condition','Subject'});
%    [pM{i},tblM{i},statsM{i},termsM{i}] = anovan(Performance,{Condition,Subject},'random',2,'varnames',{'Condition','Subject'});    
    %stu_res = assumptions_checking2(pM{i},tblM{i},statsM{i},termsM{i},Performance,Condition,Subject);
    %delete(findall(0))
    if i==6
    tableSign = multcompare(statsM{i},'CType','hsd');
    godisawoman = 1;
    end
    %% Run Friedman, cause why not
    %[pF(i), tblF{i}] = friedman(Performance_mat,1);
end