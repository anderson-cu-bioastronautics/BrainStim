% Main SR Analysis Script
% Created by Sage Sherman
% Also commented by Sage Sherman, sorry
% Last edited 10/19

%% Housekeeping and defining
clear all; close all; clc;
subID = input('What is the Subject ID? ');
TestType = input('Input the Perceptual Channel being tested: Enter \n 1 for Visual \n 2 for Auditory \n 3 for Tactile \n 4 for Vestibular \n');
SRType = input('Input the SR Type being applied: Enter \n 1 for ASR \n 2 for VSR \n');

%% Assign Variables (change code to define levels) and extract information
switch TestType
    case 1 % define information for visual test
        clusterSize=10; % for 50 trials
        if SRType == 1
            levels = [0 30:5:80];
            PlotTitle = ' ASR Visual Task'; xTitle = 'ASR Level, dB'; yTitle = 'Contrast Threshold';
        elseif SRType == 2
            levels = [0:.1:1];
            PlotTitle = ' VSR Visual Task'; xTitle = 'VSR Level, mA'; yTitle = 'Contrast Threshold';
        else
        end
    case 2 % define information for auditory task
        clusterSize=20; % for 100 trials
        if SRType == 1
            levels = [0 5.3-10:5:5.3+20 40];
            PlotTitle = ' ASR Auditory Task'; xTitle = 'ASR Level, dB'; yTitle = 'Auditory Threshold, dB';
        elseif SRType == 2
            levels = [0:.1:1];
            PlotTitle = ' VSR Auditory Task'; xTitle = 'VSR Level, mA'; yTitle = 'Auditory Threshold, dB';
        else
        end
    case 3 % define information for tactile task
        clusterSize=10;
        if SRType == 1
            levels = [0 30:5:80];
            PlotTitle = ' ASR Tactile Task'; xTitle = 'ASR Level, dB'; yTitle = 'Tactile Threshold, mm';
        elseif SRType == 2
            levels = [0:.1:1];
            PlotTitle = ' VSR Tactile Task'; xTitle = 'VSR Level, mA'; yTitle = 'Tactile Threshold, mm';
        else
        end
    case 4 % define information for vestibular task
        if SRType == 1
            levels = [0 30:5:80];
            PlotTitle = ' ASR Vestibular Task'; xTitle = 'ASR Level, dB'; yTitle = 'Vestibular Threshold, deg';
        elseif SRType == 2
            levels = [0:.1:1];
            PlotTitle = ' VSR Vestibular Task'; xTitle = 'VSR Level, mA'; yTitle = 'Vestibular Threshold, deg';
        else
        end
end

[AllData]=FindThresholds(subID,TestType,SRType,levels);

%% Plot Results
PlotData(AllData,subID,TestType,SRType,levels,PlotTitle,xTitle,yTitle,0);

%% Compare level
k=0;
while k < 1
    doStat=input('Can we compare two noise levels? (Meaning tested one level twice) \n 1 for yes \n 2 for no \n');
    switch doStat
        case 1
            Compare2Levels(subID,TestType,SRType,levels)
            k=1;
        case 2
            k=1;
    end
end

%% Further Analysis?? Bubble plots for specific noise levels
k=0;
while k < 1
    doMore=input('Do you want to look at individual plots? \n 1 for yes \n 2 for no \n');
    switch doMore
        case 1
            sinLev=input('What Noise Level do you want to observe? (Numbers only) \n');
            Int=find(levels==sinLev); % extract information of use
            if isempty(Int)==1
                while isempty(Int)==1
                    sinLev=input('ERROR, Input a proper noise level. (Numbers only) \n');
                    Int=find(levels==sinLev); % extract information of use
                end
            end
            stim=AllData{Int,1}; cor=AllData{Int,2}; lambda=AllData{Int,12};
            threshDat_bi=AllData{Int,3}; threshDat_hat=AllData{Int,4}; % does clusterSize matter? who knows?
            mu_bi=threshDat_bi(1); sigma_bi=threshDat_bi(2);
            mu_hat=threshDat_hat(1); sigma_hat=threshDat_hat(2); guessRate=0.5;
            titlestr=[num2str(subID) PlotTitle ' ' num2str(sinLev)]; xstr='Stimulus'; ystr='Probability Correct';
            plotSBC(stim,cor,lambda,mu_bi,sigma_bi,mu_hat,sigma_hat,guessRate,titlestr,xstr,ystr,clusterSize);
        case 2
            k=1;
    end
end