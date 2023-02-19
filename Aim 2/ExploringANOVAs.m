close all; clear all; clc;

%% Read in and organize data for all four trials

dataAud = readtable('DataAud.csv');
dataCAL = readtable('DataCAL.csv');
dataLPC = readtable('DataLPC.csv');
dataRMS = readtable('DataRMS.csv');
dataSmooth = readtable('DataSmooth.csv');
dataTact = readtable('DataTact.csv');
dataTimeIn = readtable('DataTimeIn.csv');
dataComp = readtable('DataComp.csv');

matAudM = []; matTactM = [];
% Average for the two perceptual readings
for i = 1:64
    datChunkA = table2array(dataAud(i*8-7:i*8,:));
    meanChunkA = [mean(datChunkA([1 5],:));mean(datChunkA([2 6],:));mean(datChunkA([3 7],:));mean(datChunkA([4 8],:))];
    datChunkT = table2array(dataTact(i*8-7:i*8,:));
    meanChunkT = [mean(datChunkT([1 5],:));mean(datChunkT([2 6],:));mean(datChunkT([3 7],:));mean(datChunkT([4 8],:))];
    if i == 1
        matAudM(1:4,:) = meanChunkA;
        matTactM(1:4,:) = meanChunkT;
    else
        matAudM(end+1:end+4,:) = meanChunkA;
        matTactM(end+1:end+4,:) = meanChunkT;
    end
end
dataAud2 = table(matAudM(:,1),matAudM(:,2),matAudM(:,3),matAudM(:,4),'VariableNames',{'Performance','Condition','Subject','Maps'});
dataTact2 = table(matTactM(:,1),matTactM(:,2),matTactM(:,3),matTactM(:,4),'VariableNames',{'Performance','Condition','Subject','Maps'});

dataAll = [dataAud2;dataCAL;dataLPC;dataRMS;dataSmooth;dataTact2;dataTimeIn];
dataFlight = [dataRMS;dataSmooth;dataTimeIn];
dataDec = [dataLPC;dataCAL];
dataPerc = [dataAud2;dataTact2];

%% Read in and organize data for all four trials

AvgdataAud = readtable('AvgDataAud.csv');
AvgdataCAL = readtable('AvgDataCAL.csv');
AvgdataLPC = readtable('AvgDataLPC.csv');
AvgdataRMS = readtable('AvgDataRMS.csv');
AvgdataSmooth = readtable('AvgDataSmooth.csv');
AvgdataTact = readtable('AvgDataTact.csv');
AvgdataTimeIn = readtable('AvgDataTimeIn.csv');
AvgdataComp = readtable('AvgDataComp.csv');

AvgdataAll = [AvgdataAud;AvgdataCAL;AvgdataLPC;AvgdataRMS;AvgdataSmooth;AvgdataTact;AvgdataTimeIn];
AvgdataFlight = [AvgdataRMS;AvgdataSmooth;AvgdataTimeIn];
AvgdataDec = [AvgdataLPC;AvgdataCAL];
AvgdataPerc = [AvgdataAud;AvgdataTact];
%% Run ANOVAs on full dataset

% % All data compiled
 [pMa,tblMa,statsMa,termsMa] = anovan(dataAll.Performance,{dataAll.Condition,dataAll.Subject,dataAll.Maps},'model',[1 0 0; 0 0 1; 1 1 0; 1 0 1],'random',2,'varnames',{'Condition','Subject','Map Difficulty'}); 
 %[pMa,tblMa,statsMa,termsMa] = anovan(dataAll.Performance,{dataAll.Condition,dataAll.Subject},'model',[1 0 ; 0 1],'random',2,'varnames',{'Condition','Subject'}); 
 %[pMa,tblMa,statsMa,termsMa] = anovan(dataAll.Performance,{dataAll.Condition,dataAll.Subject},'model',[1 0 ; 1 1],'varnames',{'Condition','Subject'}); 
% [pMa,tblMa,statsMa,termsMa] = anovan(dataAll.Performance,{dataAll.Condition,dataAll.Subject,dataAll.Maps},'model',[1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1],'varnames',{'Condition','Subject','Map Difficulty'});
 stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,dataAll.Performance,dataAll.Condition,dataAll.Subject);
% %[pMa,tblMa,statsMa,termsMa] = anovan(dataAll.Performance,{dataAll.Condition,dataAll.Subject},'model',[1 0;0 1; 1 1],'varnames',{'Condition','Subject'});
% %stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,dataAll.Performance,dataAll.Condition,dataAll.Subject);
 [AvgpMa,AvgtblMa,AvgstatsMa,AvgtermsMa] = anovan(AvgdataAll.Performance,{AvgdataAll.Condition,AvgdataAll.Subject},'model',[1 0; 1 1],'random',2,'varnames',{'Condition','Subject'});
 stu_res = assumptions_checking2(AvgpMa,AvgtblMa,AvgstatsMa,AvgtermsMa,AvgdataAll.Performance,AvgdataAll.Condition,AvgdataAll.Subject);

%% Run ANOVAs on sub datasets
for i = 1:8
    switch i
        case 1
        data = dataRMS;
        Avgdata = AvgdataRMS;       
        case 2
        data = dataTimeIn;    
        Avgdata = AvgdataTimeIn;    
        case 3
        data = dataSmooth;    
        Avgdata = AvgdataSmooth;    
        case 4
        data = dataLPC;
        Avgdata = AvgdataLPC;
        case 5
        data = dataCAL;
        Avgdata = AvgdataCAL;
        case 6
        data = dataAud;
        Avgdata = AvgdataAud;
        case 7
        data = dataTact;
        Avgdata = AvgdataTact;
        case 8 
        data = dataComp;
        Avgdata = AvgdataComp;
    end
    
%[pMa,tblMa,statsMa,termsMa] = anovan(data.Performance,{data.Condition,data.Subject,data.Maps},'model',[1 0 0; 0 0 1; 1 1 0; 1 0 1],'random',2,'varnames',{'Condition','Subject','Map Difficulty'});
 %stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,data.Performance,data.Condition,data.Subject);
 %k=1;
 %delete(findall(0))
[AvgpMa,AvgtblMa,AvgstatsMa,AvgtermsMa] = anovan(Avgdata.Performance,{Avgdata.Condition,Avgdata.Subject},'model',[1 0; 1 1],'random',2,'varnames',{'Condition','Subject'});
stu_res = assumptions_checking2(AvgpMa,AvgtblMa,AvgstatsMa,AvgtermsMa,Avgdata.Performance,Avgdata.Condition,Avgdata.Subject);
k=1;
delete(findall(0))
end

%% Run ANOVAs on subdimension dataset

 % Subdimensions
 [pMa,tblMa,statsMa,termsMa] = anovan(dataFlight.Performance,{dataFlight.Condition,dataFlight.Subject,dataFlight.Maps},'model',[1 0 0; 0 0 1; 1 1 0; 1 0 1],'random',2,'varnames',{'Condition','Subject','Map Difficulty'}); 
 stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,dataFlight.Performance,dataFlight.Condition,dataFlight.Subject);
 results = multcompare(statsMa,'Dimension',[1 3]);
 delete(findall(0))
 
 [pMa,tblMa,statsMa,termsMa] = anovan(dataDec.Performance,{dataDec.Condition,dataDec.Subject,dataDec.Maps},'model',[1 0 0; 0 0 1; 1 1 0; 1 0 1],'random',2,'varnames',{'Condition','Subject','Map Difficulty'}); 
 stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,dataDec.Performance,dataDec.Condition,dataDec.Subject);
 results = multcompare(statsMa,'Dimension',[1 3]);
 delete(findall(0))
 
 [pMa,tblMa,statsMa,termsMa] = anovan(dataPerc.Performance,{dataPerc.Condition,dataPerc.Subject,dataPerc.Maps},'model',[1 0 0; 0 0 1; 1 1 0; 1 0 1],'random',2,'varnames',{'Condition','Subject','Map Difficulty'}); 
 stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,dataPerc.Performance,dataPerc.Condition,dataPerc.Subject);
 results = multcompare(statsMa,'Dimension',[1 3]);
 delete(findall(0))
 
  % Average Subdimensions
 [pMa,tblMa,statsMa,termsMa] = anovan(AvgdataFlight.Performance,{AvgdataFlight.Condition,AvgdataFlight.Subject},'model',[1 0 ; 1 1],'random',2,'varnames',{'Condition','Subject'}); 
 stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,AvgdataFlight.Performance,AvgdataFlight.Condition,AvgdataFlight.Subject);
 [pMa,tblMa,statsMa,termsMa] = anovan(AvgdataDec.Performance,{AvgdataDec.Condition,AvgdataDec.Subject},'model',[1 0 ; 1 1],'random',2,'varnames',{'Condition','Subject'}); 
 stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,AvgdataDec.Performance,AvgdataDec.Condition,AvgdataDec.Subject);
 [pMa,tblMa,statsMa,termsMa] = anovan(AvgdataPerc.Performance,{AvgdataPerc.Condition,AvgdataPerc.Subject},'model',[1 0 ; 1 1],'random',2,'varnames',{'Condition','Subject'}); 
 stu_res = assumptions_checking2(pMa,tblMa,statsMa,termsMa,AvgdataPerc.Performance,AvgdataPerc.Condition,AvgdataPerc.Subject);
