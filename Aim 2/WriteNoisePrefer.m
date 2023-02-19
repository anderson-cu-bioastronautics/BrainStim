clear all; close all; clc;

subID = input('What is the Subject ID? ');

noise_prefer = [3 3 4 2];

    % Create a directory for the subject and their folder
currentFolder = pwd;
folder = [pwd, '\', 'Subject_Data'];
files = fullfile(folder, ['Subject_' num2str(subID)]);
    
        filename=['NoisePrefer_' num2str(subID) '.xls'];
        % Create subject's excel spreadsheet
        writeFile = fullfile(files,filename);
        xlswrite(writeFile,noise_prefer)