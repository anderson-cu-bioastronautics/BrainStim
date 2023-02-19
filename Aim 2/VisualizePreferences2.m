clear all; close all; clc;

%% Read in .mat files
%each row is a subject, each column is nGVS, AWN, MMSR, each panel is for a test (DSST, PVT, etc)

Subjects = [10:17 26:33];
numSubs = length(Subjects);
filtNum = input('How many outliers filters should occur? ');

currentFolder = pwd;
folder = [pwd, '\', 'Subject_Data'];
saveFolder = [pwd, '\', 'plots'];
mkdir(saveFolder);

for sub = 1:length(Subjects)
    files = fullfile(folder, ['Subject_' num2str(Subjects(sub))]);
    noiseData = xlsread(fullfile(files,['NoisePrefer_' num2str(Subjects(sub)) '.xls']));
    
    MP(sub) = noiseData(1); NP(sub) = noiseData(2); QP(sub) = noiseData(3);
    
    load(fullfile(files,['PerformanceData' num2str(Subjects(sub)) 'Retest.mat']));
    levelSub = SRPerfData{1,1}(7:10);
    [levSort,SortInd] = sort(levelSub); % Need to sort data in an appropriate way
    SortWant = [SortInd(1:2) SortInd(4) SortInd(3)];
    levels = levelSub(SortWant);
    
    for RetD = 2:9
        RetPerfNS(RetD-1,:) = SRPerfData{1,RetD}(7:10);
        RetPerf(RetD-1,:) = RetPerfNS(RetD-1,SortWant);
        DiffRet(RetD-1,:) = RetPerf(RetD-1,2:4)-RetPerf(RetD-1,1);
        PercRet(RetD-1,:) = DiffRet(RetD-1,:)./RetPerf(RetD-1,1)*100;
    end
    SubLevels(sub,:) = levels;
    PerfPercAll{sub,1} = PercRet;
    
    rms(sub,:) = [PercRet(1,:)];
    fracTimeIn(sub,:) = [PercRet(2,:)];
   % LZChoice(sub,:) = [PercRet(3,:)];
   % CAL(sub,:) = [PercRet(4,:)];
    Heading(sub,:) = [PercRet(3,:)];
    aud(sub,:) = [PercRet(4,:)];
    tac(sub,:) = [PercRet(5,:)];
    total(sub,:) = [PercRet(6,:)];
end

NoisePrefer = NP-QP; titCond = 'NP-QP';

%for i=1:8
for i=1:5
    base = i*numSubs-(numSubs-1);
    last = i*numSubs;
    NoisePreferAll(base:last,:) = NoisePrefer';
end


    rmsz = (rms(:,1:3)-mean(mean(rms)))./std(std(rms));
    fracTimeInz = (fracTimeIn(:,1:3)-mean(mean(fracTimeIn)))./std(std(fracTimeIn));
   % LZChoicez = (LZChoice(:,1:3)-mean(LZChoice,'all'))./std(LZChoice,0,'all');
   % CALz = (CAL(:,1:3)-mean(CAL,'all'))./std(CAL,0,'all');
    Headingz = (Heading(:,1:3)-mean(mean(Heading)))./std(std(Heading));
    audz = (aud(:,1:3)-mean(mean(aud)))./std(std(aud));
    tacz = (tac(:,1:3)-mean(mean(tac)))./std(std(tac));
  %  totalz = (total(:,1:3)-mean(total,'all'))./std(total,0,'all');


for kay = 1:3 % place all data in one array based on conditions
    %testData_mat{kay} = [rms(:,kay);fracTimeIn(:,kay);LZChoice(:,kay);CAL(:,kay);Heading(:,kay);aud(:,kay);tac(:,kay);total(:,kay)];
    %testData_mat{kay} = [rms(:,kay);fracTimeIn(:,kay);LZChoice(:,kay);CAL(:,kay);Heading(:,kay);aud(:,kay);tac(:,kay)];
    testData_mat{kay} = [rmsz(:,kay);fracTimeInz(:,kay);Headingz(:,kay);audz(:,kay);tacz(:,kay)];
    NoisePrefer_mat{kay} = NoisePreferAll;
end

for kay=1:3 % Place all data in cell arrays, maybe we'll use this for individual models
    IndtestData_cell{1,kay} = rms(:,kay);
    IndtestData_cell{2,kay} = fracTimeIn(:,kay);
    IndtestData_cell{3,kay} = Heading(:,kay);
%     IndtestData_cell{4,kay} = LZChoice(:,kay);
%     IndtestData_cell{5,kay} = CAL(:,kay);
    IndtestData_cell{4,kay} = aud(:,kay);
    IndtestData_cell{5,kay} = tac(:,kay);
 %   IndtestData_cell{6,kay} = total(:,kay);
    for lay = 1:5
        NoisePrefer_cell{lay,kay} = NoisePreferAll;
    end
end


%% Plot Data, all data including outliers
[mdlD,NoisePreferN_mat,testDataN_mat] = analyzePreferences2(NoisePrefer_mat,testData_mat);
plotPreferences2(NoisePrefer_mat,testData_mat,mdlD,numSubs,titCond,saveFolder) % Do the outlier analysis and plotting later

testDataN_mat = []; NoisePreferN_mat = []; mdlD = [];

%% Filter through data and remove outliers

for i=1:filtNum % iterate and remove outliers
    
    [mdlD,NoisePreferN_mat,testDataN_mat] = analyzePreferences2(NoisePrefer_mat,testData_mat);
    testData_mat = []; NoisePrefer_mat = []; 
    
    for kay = 1:3
        testData_mat{kay} = testDataN_mat{kay};
        NoisePrefer_mat{kay} = NoisePreferN_mat{kay};
        
        mdlDind = mdlD{kay};
        
        modelsDP(i,kay) = mdlDind.Coefficients(2,4); % Extract model significance
        mdlDind = [];
    end
    
    if i~=filtNum
        testDataN_mat = []; NoisePreferN_mat = []; mdlD = [];
    end
end

% Plot the final analyzed data
%titCond = [titCond ' NoOutliers x' num2str(filtNum)];
plotPreferences2(NoisePrefer_mat,testData_mat,mdlD,numSubs,titCond,saveFolder) % Do the outlier analysis and plotting later
