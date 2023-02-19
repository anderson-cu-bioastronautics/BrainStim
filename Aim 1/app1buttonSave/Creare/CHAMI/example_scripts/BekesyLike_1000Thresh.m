% Example script of how to use the Bekesey-Like Audiometry Exam
%
% 23 January 2013 JCW
% 29 December 2014 JCW Updated to New Specifications

%% Housekeeping:
close all; clear all; clc;

%% Inputs (this covers most user-accesable controls)
testFrequencies = [1000]; %Hz
startLevel = 40; %dB
levelUnits = 'dB SPL'; %'dB HL' or 'dB SPL'
toneType = 'sine'; %'sine', 'warble', 'noise'; %set warble and noise parameters below if you want more than the default
toneDuration = 250; %ms; duration of a presentation
toneRamp = 25; %ms; ramp-up and ramp-down duration at start and end of tone
reversalsDiscard = 2; %ignore this many reversals
reversalsKeep = 10; %test ceases after this many good reversals
incrementStart = 4; %dB, step size until "reversalsDiscard" is reached
incrementNominal = 2; %dB, setp size after "reversalsDiscard" is reached
presentationMaximum = 100; %maximum number of presentations
presentationRepetitionFrequency = 2; %Hz
maximumUnresponsivePresentations = 6; %if a subject doesn't respond at the maximum or minimum possible levels after this many presentations the exam exits

% not setting Minimum/MaximumOutputLevels because if they go beyond
% calibration min/max levels, then will error out

%% Setup the Exam:
%Start talking to the CHA:
if ~exist('chami', 'class'),
    addpath ..\src
end
cha = chami();
cha.open();

% Setup the Exam Options:
bl_exam = bekesylike_exam(); %Sets up an exam and gets the default options
bl_exam.Lstart = startLevel;
bl_exam.ToneGeneration.ToneDuration = toneDuration;
bl_exam.ToneGeneration.ToneRamp = toneRamp;
bl_exam.IncrementStart = incrementStart;
bl_exam.IncrementNominal = incrementNominal;
bl_exam.PresentationMax = presentationMaximum;
bl_exam.ReversalDiscard = reversalsDiscard;
bl_exam.ReversalKeep = reversalsKeep;
bl_exam.ToneRepetitionInterval = 1000*(1/presentationRepetitionFrequency); %the CHA wants to see ms
bl_exam.LevelUnits = levelUnits;
bl_exam.UnresponsiveMax = maximumUnresponsivePresentations;

switch toneType %sets exam parameters for the tone type specified.
    case 'sine'
        bl_exam.ToneGeneration.UseNthOctave = 0;
        bl_exam.ToneGeneration.FDevForm = 'None';
    case 'warble'
        bl_exam.ToneGeneration.UseNthOctave = 0;
        bl_exam.ToneGeneration.FDevForm = 'Sine';
        bl_exam.ToneGeneration.FDev = 5.7;
        bl_exam.ToneGeneration.FDevRate = 20;
    case 'noise'
        bl_exam.ToneGeneration.UseNthOctave = 1;
        bl_exam.ToneGeneration.FDevForm = 'None';
        bl_exam.ToneGeneration.OctaveBandSize = 3; %1/this ocataves
    otherwise
        error('Bad toneType')
end


%Prepare to loop through ears and frequencies:
thresholds = NaN.*ones(2,length(testFrequencies)); %A matrix to fill with threshold values
levels = NaN.*ones(2,length(testFrequencies),presentationMaximum ); %A matrix to fill with presentation values

%% Run the Exams:
%Loop through channels: (left / right ear):
for channel = 2 %left ear then right ear
    switch channel
        case 1
            chanName = 'Left';
            bl_exam.ToneGeneration.OutputChannel = {'HPL0' 'NONE' 'NONE' 'NONE'}; %we want the left channel First
        case 2
            chanName = 'Right';
            bl_exam.ToneGeneration.OutputChannel = {'NONE' 'HPR0' 'NONE' 'NONE'}; %we want the left channel First
    end
    
    %Loop through frequencies:
    for ff = 1:1:length(testFrequencies);
        bl_exam.F = testFrequencies(ff); %set the test frequency
        cha.queue_exam(bl_exam); %pass the exam parameters to the CHA and start the exam
        res = cha.get_exam_results(bl_exam, 10000); %Get the results; the MATLAB implmentation of this exam takes care of polling the CHA automatically
        if isnan(res.Threshold)
            warning([chanName ' Ear: ' res.ResultType ' at ' num2str(testFrequencies(ff)) ' Hz'])
        end
        thresholds(channel,ff) = res.Threshold; %Update the threshold matrix
        levels(channel,ff,1:length(res.L)) = res.L; %Update the presentation levels matrix
        bl_exam.free;  %
        pause(2)
    end
end

cha.abort_exams();
cha.close();

disp(['Test at ' num2str(thresholds(2,1)) ' dB']);

%% Plot the Results:
%Plot the threshold data:
figure;
semilogx(testFrequencies,thresholds(1,:),'bo-',testFrequencies,thresholds(2,:),'rx-','LineWidth',2)
set(gca,'yDir','reverse','FontSize',14)
ylabel('Threshold, dB HL')
xlabel('Frequency, Hz')
ylim([-10 110])

%Plot the presentation tracing:
levelsLeft = reshape(levels(1,:,:),length(testFrequencies),[]);
levelsLeft = levelsLeft';
levelsRight = reshape(levels(2,:,:),length(testFrequencies),[]);
levelsRight = levelsRight';

% figure;
% subplot(1,2,1)
% plot(levelsLeft,'LineWidth',2)
% set(gca,'FontSize',14)
% ylabel('Level, dB HL')
% xlabel('Presentation')
% legend(num2str(testFrequencies'))
% title('Left Ear')
% 
% subplot(1,2,2)
% plot(levelsRight,'LineWidth',2)
% set(gca,'FontSize',14)
% ylabel('Level, dB HL')
% xlabel('Presentation')
% title('Right Ear')
