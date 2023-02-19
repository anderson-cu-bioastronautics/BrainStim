% Main SR Analysis Script
% Created by Sage Sherman
% Also commented by Sage Sherman, sorry
% Last edited 10/19
function Main_IndSRAnalysis_post(subID,TestType,SRType,plotColor,levels,plotSave,retestAnalysis,optimal)

%% Housekeeping and defining

close all; colors;
% subID = input('What is the Subject ID? ');
% TestType = input('Input the Perceptual Channel being tested: Enter \n 1 for Visual \n 2 for Auditory \n 3 for Tactile \n 4 for Vestibular \n');
% SRType = input('Input the SR Type being applied: Enter \n 1 for ASR \n 2 for VSR \n');

% Ask about plot colors
% plotColor = input('What color do you want to make the plots?\n 0 for Black\n 1 for Red\n 2 for Blue\n 3 for Green\n 4 for Purple\n 5 for Lavender\n');


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

[InitialData]=FindThresholds(subID,TestType,SRType,levels);


%% Plot Results
PlotData(InitialData,subID,TestType,SRType,levels,PlotTitle,xTitle,yTitle,plotColor);

%% Analyze retest data
if retestAnalysis == 1
    
    % Determine how much to shift the retest plots
    if SRType == 1 % ASR
        retest_shift = 1;
        if TestType == 2
            shamRetestLevel = levels(1) - 9;
            optRetestLevel = optimal+retest_shift;
        else
            shamRetestLevel = 22;
            optRetestLevel = optimal+retest_shift;
        end
       
    elseif SRType == 2 % VSR
        retest_shift = 0.02;
        shamRetestLevel = -.08;
        optRetestLevel = optimal+retest_shift;
    else
    end
    
  
    % Generate the retest data
    RetestData = RetestAnalysis(subID,TestType,SRType,optimal);
    
    % Return to the SR curves
    figure(1);
    hold on;
    
    % Plot the retested data points
    Retest_mu = [RetestData.sham2opt2.X1,RetestData.sham2opt2.X2];
    Retest_se = [RetestData.sham2opt2.SE1,RetestData.sham2opt2.SE2];
    Retest_p = RetestData.sham2opt2.p;
    
    % Plot the retest points
    errorbar([shamRetestLevel,optRetestLevel],Retest_mu,Retest_se*1.96,'Color',colorOrder(1,:)/255);
    plot([shamRetestLevel,optRetestLevel],Retest_mu,'.','Color',colorOrder(1,:)/255,'MarkerSize',20,'LineStyle','none');

    % Annotate with the p-value
    p_str = sprintf('Retest Sham/Optimal P-value: %4.3f',Retest_p);
    boxPos = [0.2 0.2 0.1 0.6];
    annotation('textbox',boxPos,'String',p_str,'FitBoxToText','on','EdgeColor','none','BackgroundColor','none');
end

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

%% Save the Plots and data

% Do you want to save the plots? 1 for yes, 0 for no
% plotSave = input('Do you want to save the plots? Enter \n 1 for Yes \n 0 for No\n');

% The directory to save the plots
saveFigDir = ['RA1Data/',num2str(subID),'/Figures/'];

% The data file name
switch TestType
    case 1 % Visual
        TaskName = 'Visual_';
    case 2 % Auditory
        TaskName = 'Auditory_';
    case 3 % Tactile
        TaskName = 'Tactile_';
    case 4 %  Vestibular
        TaskName = 'Vestibular_';
end

if SRType ==1
    SRName = 'ASR';
elseif SRType ==2
    SRName = 'VSR';
end

file2save = ['RA1Data/',num2str(subID),'/Subject_',num2str(subID),'_',TaskName,SRName,'_Initial'];

% Save the data
save(file2save,'InitialData');

if retestAnalysis == 1
    file2save = ['RA1Data/',num2str(subID),'/Subject_',num2str(subID),'_',TaskName,SRName,'_Retest'];
    save(file2save,'RetestData');
end

% Save the plots, if that's what you want to do
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
    figName = ['Subject ' num2str(subID) PlotTitle ' White Background'];
    saveStrFIG = strcat(saveFigDir,figName);
    saveStrPNG = strcat(saveFigDir,figName);
    saveStrPDF = strcat(saveFigDir,figName);
    saveas(gcf,[saveStrFIG '.fig']);
    print(saveStrPNG,'-dpng','-r100');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 4.5]);
    set(gcf,'PaperOrientation','landscape');
    print(saveStrPDF,'-dpdf','-r600');
    
%     % Make the figure with the poster background
%     set(gca,'Color',[250 245 252]/255);
%     set(gcf,'Color',[250 245 252]/255);
%     
%     % Save the figure with the default white background
%     figName = ['Subject ' num2str(subID) PlotTitle ' Poster Background'];
%     saveStrFIG = strcat(saveFigDir,figName);
%     saveStrPNG = strcat(saveFigDir,figName);
%     saveStrPDF = strcat(saveFigDir,figName);
%     saveas(gcf,[saveStrFIG '.fig']);
%     print(saveStrPNG,'-dpng','-r100');
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 4.5]);
%     set(gcf,'PaperOrientation','landscape');
%     print(saveStrPDF,'-dpdf','-r600');
    
%     % Do the same thing for figure 2
%     figure(2);
%     
%     % Make the figure pretty, with white background
%     set(gca,'FontName','Arial');
%     set(gca,'Color',[1 1 1]);
%     set(gcf,'Color',[1 1 1]);
%     set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 6 4]);
%     
%     % Save the figure with the default white background
%     figName = ['Subject ' num2str(subID) PlotTitle ' White Background'];
%     saveStrFIG = strcat(saveFigDir,'FIG/',figName);
%     saveStrPNG = strcat(saveFigDir,'PNG/',figName);
%     saveStrPDF = strcat(saveFigDir,'PDF/',figName);
%     saveas(gcf,[saveStrFIG '.fig']);
%     print(saveStrPNG,'-dpng','-r100');
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 4.5]);
%     set(gcf,'PaperOrientation','landscape');
%     print(saveStrPDF,'-dpdf','-r600');
%     
%     % Make the figure with the poster background
%     set(gca,'Color',[250 245 252]/255);
%     set(gcf,'Color',[250 245 252]/255);
%     
%     % Save the figure with the default white background
%     figName = ['Subject ' num2str(subID) PlotTitle ' Poster Background'];
%     saveStrFIG = strcat(saveFigDir,'FIG/',figName);
%     saveStrPNG = strcat(saveFigDir,'PNG/',figName);
%     saveStrPDF = strcat(saveFigDir,'PDF/',figName);
%     saveas(gcf,[saveStrFIG '.fig']);
%     print(saveStrPNG,'-dpng','-r100');
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 4.5]);
%     set(gcf,'PaperOrientation','landscape');
%     print(saveStrPDF,'-dpdf','-r600');
end