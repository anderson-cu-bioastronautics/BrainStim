function [writeFile]=OpenWriteFile(subNum,filename)

% Create a directory for the subject
currentFolder = pwd;
SubjectFolder = fullfile(currentFolder,'RA1Data',num2str(subNum));
mkdir(SubjectFolder);

% Create subject's excel spreadsheet
writeFile = fullfile(SubjectFolder,filename);
Trial=0;Stim=0;Response=0;Correct=0;
headers=table(Trial,Stim,Response,Correct);
writetable(headers,writeFile);
