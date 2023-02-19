%% Analyze all the subjects
clc;
clear;
close all;
tic;

% Ask about plot colors
plotColor = input('What color do you want to make the plots?\n 0 for Black\n 1 for Red\n 2 for Blue\n 3 for Green\n 4 for Purple\n 5 for Lavender\n');

% Save the plots!
plotSave = true;

% What subjects are we analyzing?
subjects = [1 2 7]; %????

% Tasks [Visual, Auditory, Tactile]
tasks = [1 2 3];

% SR types [ASR, VSR]
SRTypes = [1 2];

% VSR Levels
VSRLevels = [0:0.1:1];

% ASR Levels
ASRLevels = [0 30:5:80];

% Subject [] Auditory ASR Levels
sub1AudASRLevels = [0 3 8 13 18 23 28 40];

% Subject [] Auditory ASR Levels
sub2AudASRLevels = [0 3 8 13 23 28 33 40];

% Subject [] Auditory ASR Levels
sub7AudASRLevels = [-4 0 1 6 11 16 21 26 40];


%% Loop through and do the thing!

% First, categorize by SR type
for SRType = 1:length(SRTypes)
    
    % Second, categorize by task type
    for TestType = 1:length(tasks)
        
        % Third, categorize by subject
        for subID = subjects
            
            % Determine the SR levels
            
            % Are we working with ASR?
            if SRType == 1
                
                % Is this the auditory task?
                if TestType ==2
                    
                    % The levels are those applied uniquely to each
                    % subject!
                    if subID == subjects(1)
                        levels = sub1AudASRLevels;
                    elseif subID == subjects(2)
                        levels = sub2AudASRLevels;
                    else
                        levels =sub7AudASRLevels;
                    end
                else
                    % If not auditory task, use the default ASR levels.
                    levels = ASRLevels;
                end
            else
                % If VSR, use default VSR levels.
                levels = VSRLevels;
            end
            
            % Now, do the analysis, create the figures, and save them!
            Main_IndSRAnalysis_new(subID,TestType,SRType,plotColor,levels,plotSave);
            
        end
    end
end

programRunTime = toc;
fprintf('Program Run Time: %4.1f\n',programRunTime);
            