%% Script to create and store randomized list

clear all; close all; clc;
Subject = input('What is the Subject ID? ');

% Subject file information
if Subject >10
    Subject_str = num2str(Subject);
else
    Subject_str = ['0',num2str(Subject)];
end 
%% Determine and identify the task for naming convention

% Same ASR and VSR Noise Levels for all tasks, so that's nice
levels = [0.2:0.3:0.8 40:15:70];

rand_noise=randperm(length(levels));
noise_order=levels(rand_noise)

% Create a directory for the subject and their folder
currentFolder = pwd;
SubjectFolder = fullfile(currentFolder,'Subject_Data',['Subject_' Subject_str]);
mkdir(SubjectFolder);

filename=['TestingOrder_' Subject_str '.xls'];

% Create subject's excel spreadsheet
writeFile = fullfile(SubjectFolder,filename);
xlswrite(writeFile,noise_order)

% %% Randomize the maps
% 
% % Define Map vectors
% easyO = [10 16 19 23 27 31 36];
% medO = [11 12 14 15 18 21 24 25 26 28 30 32 34 37];
% hardO = [13 17 20 22 29 33 35];
% 
% easyR = [38 42 48 53];
% medR = [40 41 43 44 46 49 50 51];
% hardR = [39 45 47 52];
% 
% %randomize map vectors
% easyO_randvect = randperm(length(easyO));
% medO_randvect = randperm(length(medO));
% hardO_randvect = randperm(length(hardO));
% 
% easyR_randvect = randperm(length(easyR));
% medR_randvect = randperm(length(medR));
% hardR_randvect = randperm(length(hardR));
% 
% Orig_mat = [easyO(easyO_randvect); medO(medO_randvect(1:7)); medO(medO_randvect(8:14)); hardO(hardO_randvect)];
% Ret_mat = [easyR(easyR_randvect); medR(medR_randvect(1:4)); medR(medR_randvect(5:8)); hardR(hardR_randvect)];
% 
% % randomize easy - hard presentation of the map matrices
% for i=1:7
%    randCol = randperm(4);
%    randO_mat(:,i) = Orig_mat(randCol,i);
% end
% 
% for i =1:4
%    randCol = randperm(4);
%    randR_mat(:,i) = Ret_mat(randCol,i);    
% end
% 
% map_order = [randO_mat randR_mat];
% 
% filenameM=['MapOrder_' num2str(subID) '.xls'];
% 
% % Create subject's excel spreadsheet
% writeFile = fullfile(SubjectFolder,filenameM);
% xlswrite(writeFile,map_order)