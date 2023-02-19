%% Script to create and store randomized list

clear all; close all; clc;
subID = input('What is the Subject ID? ');
TestType = input('Input the Perceptual Channel being tested: Enter \n 1 for Visual \n 2 for Auditory \n 3 for Tactile \n 4 for Vestibular \n');
SRType = input('Input the SR Type being applied: Enter \n 1 for ASR \n 2 for VSR \n');

if SRType == 1
    SRname='ASR';
elseif SRType == 2
    SRname='VSR';    
end
%% Determine and randomize the name of certain levels
switch TestType
    case 1 % define information for visual test
        Test='Visual';
        if SRType == 1
            levels = [0 30:5:80];
        elseif SRType == 2
            levels = [0:.1:1];
        else
        end
    case 2 % define information for auditory task
        Test='Auditory';
        if SRType == 1
            levels = [0 -3:5:27 40];
        elseif SRType == 2
            levels = [0:.1:1];
        else
        end
    case 3 % define information for tactile task
        Test='Tactile';
        if SRType == 1
            levels = [0 30:5:80];
        elseif SRType == 2
            levels = [0:.1:1];
        else
        end
    case 4 % define information for vestibular task
        Test='Vestibular';
        if SRType == 1
            levels = [0 30:5:80];
        elseif SRType == 2
            levels = [0:.1:1];
        else
        end
end

rand_noise=randperm(length(levels));
noise_order=levels(rand_noise)

% Create a directory for the subject
currentFolder = pwd;
SubjectFolder = fullfile(currentFolder,'RA1Data',num2str(subID));
mkdir(SubjectFolder);
filename=['TestingOrder_' SRname '_' num2str(subID) '_' Test '.xls'];

% Create subject's excel spreadsheet
writeFile = fullfile(SubjectFolder,filename);
xlswrite(writeFile,noise_order)