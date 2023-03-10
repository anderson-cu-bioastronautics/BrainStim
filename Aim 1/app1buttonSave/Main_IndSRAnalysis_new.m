% Main SR Analysis Script
% Created by Sage Sherman
% Also commented by Sage Sherman, sorry
% Last edited 9/10/20 by Abby Durell


function Main_IndSRAnalysis_new(subID,TestType,SRType,plotColor,levels,plotSave)
clear global; close all; clc; colors;
tic

global TestType SRType levels thisSubjectMeasuredThresholds method widthThreshold1 widthThreshold2

method = 'Derivative';

%% Housekeeping and defining


% subID = input('What is the Subject ID? ');
% TestType = input('Input the Perceptual Channel being tested: Enter \n 1 for Visual \n 2 for Auditory \n 3 for Tactile \n 4 for Vestibular \n');
% SRType = input('Input the SR Type being applied: Enter \n 1 for ASR \n 2 for VSR \n');

% Ask about plot colors
% plotColor = input('What color do you want to make the plots?\n 0 for Black\n 1 for Red\n 2 for Blue\n 3 for Green\n 4 for Purple\n 5 for Lavender\n');
 
% ASR Auditory has masking
maskingExpected = (TestType == 2 && SRType == 1);

%% Assign Variables (change code to define levels) and extract information
switch TestType
    case 1 % define information for visual test
        clusterSize=10; % for 50 trials
        if SRType == 1
%             levels = [0 30:5:80];
            PlotTitle = ' ASR Visual Task'; xTitle = 'ASR Level, dB'; yTitle = 'Contrast Threshold';
        elseif SRType == 2
%             levels = [0:.1:1];
            PlotTitle = ' VSR Visual Task'; xTitle = 'VSR Level, mA'; yTitle = 'Contrast Threshold';
        else
        end
    case 2 % define information for auditory task
        clusterSize=20; % for 100 trials
        if SRType == 1
%             levels = input('Input the vector of Auditory Levels Tested:\n');
            PlotTitle = ' ASR Auditory Task'; xTitle = 'ASR Level, dB'; yTitle = 'Auditory Threshold, dB SPL';
        elseif SRType == 2
%             levels = [0:.1:1];
            PlotTitle = ' VSR Auditory Task'; xTitle = 'VSR Level, mA'; yTitle = 'Auditory Threshold, dB SPL';
        else
        end
    case 3 % define information for tactile task
        clusterSize=10;
        if SRType == 1
%             levels = [0 30:5:80];
            PlotTitle = ' ASR Tactile Task'; xTitle = 'ASR Level, dB'; yTitle = 'Tactile Threshold, Displacement Fraction';
        elseif SRType == 2
%             levels = [0:.1:1];
            PlotTitle = ' VSR Tactile Task'; xTitle = 'VSR Level, mA'; yTitle = 'Tactile Threshold, Displacement Fraction';
        else
        end
    case 4 % define information for vestibular task
        if SRType == 1
%             levels = [0 30:5:80];
            PlotTitle = ' ASR Vestibular Task'; xTitle = 'ASR Level, dB'; yTitle = 'Vestibular Threshold, deg';
        elseif SRType == 2
%             levels = [0:.1:1];
            PlotTitle = ' VSR Vestibular Task'; xTitle = 'VSR Level, mA'; yTitle = 'Vestibular Threshold, deg';
        else
        end
       
end

%% Find Thresholds
[AllData,thisSubjectMeasuredThresholds]=FindThresholds(subID,TestType,SRType,levels);

%% Curve Fit Setup
[allUnderlyingThresholds,allMeasuredThresholds,A0Guess,lambdaGuess,w0Guess,mGuess,sGuess,fGuess,BGuess,B_sd,Ymin,Ymax,underlyingCurveParameters ,UpperBounds,LowerBounds,options,...
    guess6,BTL,MST] = curveFitSetup(SRType,TestType,maskingExpected,levels,thisSubjectMeasuredThresholds);

%% Curve Fitting Implementation
 
underlyingShamThreshold = MST;  
 
[curveFitParameters,A0_fit,lambda_fit,w0_fit,m_fit,s_fit,f_fit,B_fit,X,r,SRDip,mask,SRCurveFit]...
    = Curve_fitting_implementation_and_guesses(underlyingCurveParameters,maskingExpected,levels,thisSubjectMeasuredThresholds,LowerBounds,UpperBounds,...
    options,underlyingShamThreshold,A0Guess,lambdaGuess,w0Guess,mGuess,sGuess,fGuess,BGuess,B_sd,guess6,BTL,MST);

%% Determine the retest noise levels
    
    % Find the curve fit minimum
    curveFitMinimum = min(X(find(SRCurveFit == min(SRCurveFit))));
  
    if SRType == 1 && TestType == 2
        % Find the rounded
    roundedCFIndex = 1+find(abs(levels(2:end-1) - curveFitMinimum) == min(abs(levels(2:end-1) - curveFitMinimum)));
    
    % Find the closest initial swath level
    roundedCFMinimum = levels(roundedCFIndex);
    
    % Find the retest indices
    if roundedCFIndex == 2
        
        % The retest indices are the best and 2 above
        retestIndices = [1,2,3,4];
        
        % Single retest indices
        singleRetestIndices = [3,4];
        
        doubleRetestIndices = [1,2];
        
        bestRetestIndex = 2;
        rawDoubleIndices = [1,2];
        rawSingleIndices = [3,4];
        
    elseif roundedCFIndex == length(levels)-1
        
        % The retest indices are the maximum SR noise level and 2 below
        retestIndices = [1, roundedCFIndex + [-2, -1, 0]];
        
        % Single retest indices
        singleRetestIndices = roundedCFIndex + [-2, -1];
        
        bestRetestIndex = 4;
        
        rawDoubleIndices = [1,4];
        rawSingleIndices = [2,3];
        
    else
        
        % The retest indices are the rounded CF minimum plus and minus one
        retestIndices = [1, roundedCFIndex + [-1, 0, 1]];
        
        % Single retest indices
        singleRetestIndices = roundedCFIndex + [-1, 1];
            
        bestRetestIndex = 3;
        
        rawDoubleIndices = [1,3];
        rawSingleIndices = [2,4];
    end
    
    retestLevels = levels(retestIndices);
    

    else
    % Find the rounded
    roundedCFIndex = 1+find(abs(levels(2:end) - curveFitMinimum) == min(abs(levels(2:end) - curveFitMinimum)));
    
    % Find the closest initial swath level
    roundedCFMinimum = levels(roundedCFIndex);
    
    % Find the retest indices
    if roundedCFIndex == 2
        
        % The retest indices are the best and 2 above
        retestIndices = [1,2,3,4];
        
        % Single retest indices
        singleRetestIndices = [3,4];
        
        doubleRetestIndices = [1,2];
        
        bestRetestIndex = 2;
        rawDoubleIndices = [1,2];
        rawSingleIndices = [3,4];
        
    elseif roundedCFIndex == length(levels)
        
        % The retest indices are the maximum SR noise level and 2 below
        retestIndices = [1, roundedCFIndex + [-2, -1, 0]];
        
        % Single retest indices
        singleRetestIndices = roundedCFIndex + [-2, -1];
        
        bestRetestIndex = 4;
        
        rawDoubleIndices = [1,4];
        rawSingleIndices = [2,3];
        
    else
        
        % The retest indices are the rounded CF minimum plus and minus one
        retestIndices = [1, roundedCFIndex + [-1, 0, 1]];
        
        % Single retest indices
        singleRetestIndices = roundedCFIndex + [-1, 1];
            
        bestRetestIndex = 3;
        
        rawDoubleIndices = [1,3];
        rawSingleIndices = [2,4];
    end
    end
    retestLevels = levels(retestIndices);
    
     rand_noise_retest=retestLevels(randperm(length(retestLevels)))
%     
     %noise_order_retest=retestLevels(rand_noise_retest)
    %fprintf('The Retest Values are %.2f, %.2f, %.2f, %.2f \n',retestLevels(1),retestLevels(2),retestLevels(3),retestLevels(4))
    
    
    %% Make sure to randomize the retest levels and say how many trials to do at each!
    
    % Visual/tactile tasks: 100 at sham/best, 50 at neighbors
    % Auditory task: 175 at sham/best, 88 at neighbors
    % ASR Auditory: TBD, so use auditory task for now
%% Plot Results
PlotData_new(AllData,subID,TestType,SRType,PlotTitle,xTitle,yTitle,plotColor,X,levels,thisSubjectMeasuredThresholds,SRCurveFit,...
    Ymin,Ymax);

%% Compare level
% k=0;
% while k < 1
%     doStat=input('Can we compare two noise levels? (Meaning tested one level twice) \n 1 for yes \n 2 for no \n');
%     switch doStat
%         case 1
%             Compare2Levels(subID,TestType,SRType,levels)
%             k=1;
%         case 2
%             k=1;
%     end
% end

%% Further Analysis?? Bubble plots for specific noise levels
% k=0;
% while k < 1
%     doMore=input('Do you want to look at individual plots? \n 1 for yes \n 2 for no \n');
%     switch doMore
%         case 1
%             sinLev=input('What Noise Level do you want to observe? (Numbers only) \n');
%             Int=find(levels==sinLev); % extract information of use
%             if isempty(Int)==1
%                 while isempty(Int)==1
%                     sinLev=input('ERROR, Input a proper noise level. (Numbers only) \n');
%                     Int=find(levels==sinLev); % extract information of use
%                 end
%             end
%             stim=AllData{Int,1}; cor=AllData{Int,2}; lambda=AllData{Int,12};
%             threshDat_bi=AllData{Int,3}; threshDat_hat=AllData{Int,4}; % does clusterSize matter? who knows?
%             mu_bi=threshDat_bi(1); sigma_bi=threshDat_bi(2);
%             mu_hat=threshDat_hat(1); sigma_hat=threshDat_hat(2); guessRate=0.5;
%             titlestr=[num2str(subID) PlotTitle ' ' num2str(sinLev)]; xstr='Stimulus'; ystr='Probability Correct';
%             plotSBC(stim,cor,lambda,mu_bi,sigma_bi,mu_hat,sigma_hat,guessRate,titlestr,xstr,ystr,clusterSize);
%         case 2
%             k=1;
%     end
% end

%% Save the Plots

% Do you want to save the plots? 1 for yes, 0 for no
% plotSave = input('Do you want to save the plots? Enter \n 1 for Yes \n 0 for No\n');

% The directory to save the plots
saveDir = ['RA1Data/Figures/'];


if plotSave == true
    
    % Start by assuming that this directory is not where we want to save the
    % plots
%     changeDir = 1;
%     while changeDir ==1
%         % Display the directory where the plots will save
%         fprintf('The current save directory is:\n  ');
%         disp(saveDir);
%         
%         % Ask if the user wants to change the save directory
%         changeDir = input('\nDo you want to change the save directory?\n 1 for Yes\n 0 for No');
%         
%         % If they want to change the directory, ask for the new directory
%         if changeDir ==1
%             saveDir = input('Where would you like to save the plots?','s');
%         
%         % Otherwise, the directory does not need to be changed.
%         else
%             changeDir = 0;
%         end
%     end
    
    
    % Make figure 1 the working figure
    figure(1);
    
    % Make the figure pretty, with white background
    set(gca,'FontName','Arial');
    set(gca,'Color',[1 1 1]);
    set(gcf,'Color',[1 1 1]);
    set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 6 4]);
    
    % Save the figure with the default white background
    figName = ['Subject ' num2str(subID) PlotTitle ' SR Curve White Background'];
    saveStrFIG = strcat(saveDir,'FIG/SRCurves/White/',figName,'.fig');
    saveStrPNG = strcat(saveDir,'PNG/SRCurves/White/',figName);
    saveStrPDF = strcat(saveDir,'PDF/SRCurves/White/',figName);
    saveas(gcf,saveStrFIG);
    print(saveStrPNG,'-dpng','-r100');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 4.5]);
    print(saveStrPDF,'-dpdf','-r600');
    
    % Make the figure with the poster background
    set(gca,'Color',[250 245 252]/255);
    set(gcf,'Color',[250 245 252]/255);
    
    % Save the figure with the poster background
    figName = ['Subject ' num2str(subID) PlotTitle ' SR Curve Poster Background'];
    saveStrFIG = strcat(saveDir,'FIG/SRCurves/Poster/',figName,'.fig');
    saveStrPNG = strcat(saveDir,'PNG/SRCurves/Poster',figName);
    saveStrPDF = strcat(saveDir,'PDF/SRCurves/Poster',figName);
    saveas(gcf,[saveStrFIG '.fig']);
    print(saveStrPNG,'-dpng','-r100');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 4.5]);
    print(saveStrPDF,'-dpdf','-r600');
    
    % Do the same thing for figure 2
    figure(2);
    
    % Make the figure pretty, with white background
    set(gca,'FontName','Arial');
    subplot(1,2,1);
    set(gca,'Color',[1 1 1]);
    subplot(1,2,2);
    set(gca,'Color',[1 1 1]);
    set(gcf,'Color',[1 1 1]);
    set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 6 4]);
    
    % Save the figure with the default white background
    figName = ['Subject ' num2str(subID) PlotTitle ' Psychometrics White Background'];
    saveStrFIG = strcat(saveDir,'FIG/PsychometricCurves/White/',figName);
    saveStrPNG = strcat(saveDir,'PNG/PsychometricCurves/White/',figName);
    saveStrPDF = strcat(saveDir,'PDF/PsychometricCurves/White/',figName);
    saveas(gcf,[saveStrFIG '.fig']);
    print(saveStrPNG,'-dpng','-r100');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 4.5]);
    print(saveStrPDF,'-dpdf','-r600');
    
    % Make the figure with the poster background
    subplot(1,2,1);
    set(gca,'Color',[250 245 252]/255);
    subplot(1,2,2);
    set(gca,'Color',[250 245 252]/255);
    set(gca,'Color',[250 245 252]/255);
    set(gcf,'Color',[250 245 252]/255);
    
    % Save the figure with the default white background
    figName = ['Subject ' num2str(subID) PlotTitle ' Psychometrics Poster Background'];
    saveStrFIG = strcat(saveDir,'FIG/PsychometricCurves/Poster/',figName);
    saveStrPNG = strcat(saveDir,'PNG/PsychometricCurves/Poster/',figName);
    saveStrPDF = strcat(saveDir,'PDF/PsychometricCurves/Poster/',figName);
    saveas(gcf,[saveStrFIG '.fig']);
    print(saveStrPNG,'-dpng','-r100');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 4.5]);
    print(saveStrPDF,'-dpdf','-r600');
end

 RunTime = toc;
 RunTime = RunTime/60;
 fprintf('RunTime= %.2f minutes',RunTime)