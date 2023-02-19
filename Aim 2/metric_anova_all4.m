%% ANOVA for Metrics Across All Maps
% Y. Shen
% 18 August 2021

%% Load data
clear all; clc; %delete(findall(0));
SubjectsNAA = 10:17;
SubjectsAA = [26:29 31:33];
Subjects = [SubjectsNAA SubjectsAA];
sizeOrg = 6;
metrics = {'RMS', 'Stick Input', 'Smoothness', 'LP Choice', 'Crash Abort Land', 'Auditory Reaction', 'Tactile Reaction', 'Total Score', 'CAL Score', 'PreChoice', 'PostChoice'};
N_metrics = 11;
conds = {'sham' 'nGVS' 'AWN' 'MMSR'};
N_conds = 4;
countNAA = 1;
countAA = 1;

count = 1;
MapMat_unsorted = [1 3 2 2; 1 2 2 3; 2 3 1 2; 2 2 3 1];

for Subject = [SubjectsNAA SubjectsAA]
    
    currentFolder = pwd;
    
    % Subject file information
    if Subject >9
        Subject_str = num2str(Subject);
    else
        Subject_str = ['0',num2str(Subject)];
    end
    
    
    load(['Subject_Data/Subject_' Subject_str '/PerformanceData' Subject_str 'Retest.mat']);
    load(['Subject_Data/Subject_' Subject_str '/RankedData' Subject_str 'Retest.mat']);
    %load(['RankedData' Subject_str '.mat']);
    
    Performance_Data1 = SRPerfData;
    Performance_Data = SRRankData;
    
    SRlevels = Performance_Data{1,1};
    rms_p = 1-Performance_Data{1,2}./16;
    rms_p(isnan(rms_p))=0;
    TimeIn_p = 1-Performance_Data{1,3}./16;
    TimeIn_p(isnan(TimeIn_p))=0;
    Choice_P_p = Performance_Data{1,4}./16;
    Choice_P_p(isnan(Choice_P_p))=1;
    CAL_p = Performance_Data{1,5}./16;
    CAL_p(isnan(CAL_p))=1;
    Head_Change_p = 1-Performance_Data{1,6}./16;
    Head_Change_p(isnan(Head_Change_p))=0;
    audStop_p = 1-Performance_Data{1,7}./32;
    audStop_p(isnan(audStop_p))=0;
    tacStop_p = 1-Performance_Data{1,8}./32;
    tacStop_p(isnan(tacStop_p))=0;
    Total_p = 1/3*(rms_p + TimeIn_p + Head_Change_p) + 1/2*(Choice_P_p + CAL_p) + 1/2*(audStop_p(1:4,:) + tacStop_p(1:4,:) + audStop_p(5:8,:) + tacStop_p(5:8,:));
    PerfData = {SRlevels; rms_p; TimeIn_p;  Head_Change_p; Choice_P_p; CAL_p; audStop_p; tacStop_p; Total_p};
    
    Order = SRlevels;
    Order2 = Order((end-3):end); %retest
    
    RetData = PerfData;%(:, (sizeOrg+1):end);
    
    %   OrigData = PerfData(:, 1:sizeOrg);
    
    [RetOrder, sorted] = sort(Order2, 'ascend'); %will sort to sham, vsr, mmsr, asr
    MapMat_sorted = MapMat_unsorted(:, sorted);
    MapMat_sorted(:,3:4) = flip(MapMat_sorted(:,3:4),2);
    for q = 1:9
        RetDataSorted{count,q} = RetData{q}(:, sorted)-mean(RetData{q}(:));
        RetDataSorted{count,q}(:,3:4) = flip(RetDataSorted{count,q}(:,3:4),2);
    end
    SubDataMat{count} = ones(4)*Subject;
    CondDataMat{count} = [ones(4,1)*1 ones(4,1)*2 ones(4,1)*3 ones(4,1)*4];
    mapDataMat{count} = MapMat_sorted;
    SubDataMatPer{count} = ones(8,4)*Subject;
    CondDataMatPer{count} = [ones(8,1)*1 ones(8,1)*2 ones(8,1)*3 ones(8,1)*4];
    mapDataMatPer{count} = [MapMat_sorted; MapMat_sorted];
    count = count+1;
    %PercentData = []; RetDataSorted = []; PerfData = []; Performance_Data =[];
end

for q = 2:9
    AllPerfMat{q-1} = cat(1,RetDataSorted{:,q});
    if q == 7 || q == 8
        AllPerfarr{q-1} = reshape(AllPerfMat{q-1},[128*4,1]);
    else
        AllPerfarr{q-1} = reshape(AllPerfMat{q-1},[64*4,1]);
    end
end
AllSubMat = cat(1,SubDataMat{:});
AllCondMat = cat(1,CondDataMat{:});
AllMapMat = cat(1,mapDataMat{:});
AllMapMatPer = cat(1,mapDataMatPer{:});
AllSubMatPer = cat(1,SubDataMatPer{:});
AllCondMatPer = cat(1,CondDataMatPer{:});
AllSubarr = reshape(AllSubMat,[64*4,1]);
AllCondarr = reshape(AllCondMat,[64*4,1]);
AllMaparr = reshape(AllMapMat,[64*4,1]);
AllMaparrPer = reshape(AllMapMatPer,[128*4,1]);
AllSubarrPer = reshape(AllSubMatPer,[128*4,1]);
AllCondarrPer = reshape(AllCondMatPer,[128*4,1]);

%% Check ANOVA Assumptions

% Test for normality

fprintf('AllSubjects:\n');
for i = 1:N_metrics
    [h, p] = kstest(AllPerfarr{i});
    fprintf('%s - H = %g, p = %g\n', metrics{i}, h, p);
end


%% ANOVAs

fprintf('\nRMANOVAs\n');
fprintf('==========\n');
rm = cell(N_metrics, 1);
for i = 1:N_metrics
    Performance = AllPerfarr{i};
    PerfMat = AllPerfMat{i};
    
    if i == 6 || i==7
    Condition = AllCondarrPer;
    Subject = AllSubarrPer;
    Maps = AllMaparrPer;
    else
    Condition = AllCondarr;
    Subject = AllSubarr;
    Maps = AllMaparr;
    end
    %[pA,tblA,statsA,termsA] = anovan(Performance,{Condition,Age,Subject},'model',[1,0,0; 0,1,0; 0,0,1],'nested',[0,0,0; 0,0,0; 0,1,0],'random',3,'varnames',{'Condition','Age','Subject'});
%     if i==8
%         [pA,tblA,statsA,termsA] = anovan(Performance,{Condition,Age,Subject},'model',[1,0,0; 0,1,0; 0,0,1],'nested',[0,0,0; 0,0,0; 0,1,0],'random',3,'varnames',{'Condition','Age','Subject'});
%     end
    %Condition = categorical(Condition);
    %Subject = categorical(Subject);
    %Maps = categorical(Maps);
%metrics = {'RMS', 'Stick Input', 'Smoothness', 'LP Choice', 'Crash Abort Land', 'Auditory Reaction', 'Tactile Reaction', 'Total Score'};
      
    tblDat = table(Performance,Condition,Subject,Maps);
    if i == 1
    writetable(tblDat,'DataRMS.csv')
    elseif i == 2
    writetable(tblDat,'DataTimeIn.csv')
    elseif i == 3
    writetable(tblDat,'DataSmooth.csv')
    elseif i == 4
    writetable(tblDat,'DataLPC.csv')
    elseif i == 5
    writetable(tblDat,'DataCAL.csv')
    elseif i == 6
    writetable(tblDat,'DataAud.csv')
    elseif i == 7
    writetable(tblDat,'DataTact.csv')
    elseif i == 8
    writetable(tblDat,'DataComp.csv')
    elseif i == 9
    writetable(tblDat,'DataCALScore.csv')
    elseif i == 10
    writetable(tblDat,'preChoice.csv')
    elseif i == 11
    writetable(tblDat,'postChoice.csv')
    end
    
    [pMa{i},tblMa{i},statsMa{i},termsMa{i}] = anovan(Performance,{Condition,Subject,Maps},'model',[1 0 0; 0 1 0; 0 0 1; 1 1 0],'random',2,'varnames',{'Condition','Subject','Maps'});
   % [pMa{i},tblMa{i},statsMa{i},termsMa{i}] = anovan(Performance,{Condition,Subject},'model',[1 0; 0 1 ; 1 1 ],'random',2,'varnames',{'Condition','Subject'});
   %[pMa{i},tblMa{i},statsMa{i},termsMa{i}] = anovan(Performance,{Condition,Subject},'random',2,'varnames',{'Condition','Subject'});
   % [pMA{i},tblMA{i},statsMA{i}] = anova2(PerfMat,4)
   % [pM{i},tblM{i},statsM{i}] = friedman(PerfMat,4);
    hmm=1;
    %[pM{i},tblM{i},statsM{i},termsM{i}] = anovan(Performance,{Condition,Subject},'random',2,'varnames',{'Condition','Subject'});
    stu_res = assumptions_checking2(pMa{i},tblMa{i},statsMa{i},termsMa{i},Performance,Condition,Subject);
    hmm = 1;
    %delete(findall(0))
%     if i==6
%         tableSign = multcompare(statsM{i},'CType','hsd');
%         godisawoman = 1;
%     end
    %% Run Friedman, cause why not
    %[pF(i), tblF{i}] = friedman(Performance_mat,1);
end