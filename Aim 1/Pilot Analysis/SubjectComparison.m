% Function to analyze SR thresholds across all subject data

clear all; close all; clc;
%subID = input('What is the Subject ID? ');
TestType = input('Input the Perceptual Channel being tested: Enter \n 1 for Visual \n 2 for Auditory \n 3 for Tactile \n 4 for Vestibular \n');
SRType = input('Input the SR Type being applied: Enter \n 1 for ASR \n 2 for VSR \n');
colors=['k' 'b' 'r' 'm' 'c'];

%% Assign Variables (change code to define levels) and extract information
switch TestType
    case 1 % define information for visual test
        if SRType == 1
            levels = [0 30:5:80];
            PlotTitle = ' ASR Visual Task'; xTitle = 'ASR Level, dB'; yTitle = 'Contrast Threshold';
            subjects = [1002 1003 1006 1007 1008];
        elseif SRType == 2
            levels = [0:.1:1];
            PlotTitle = ' VSR Visual Task'; xTitle = 'VSR Level, mA'; yTitle = 'Contrast Threshold';
            subjects = [1002 1003 1004];
        else
        end
    case 2 % define information for auditory task
        if SRType == 1
            levels = [-5 0 0.1 5:5:40];
            PlotTitle = ' ASR Auditory Task'; xTitle = 'ASR Level, dB'; yTitle = 'Auditory Threshold, dB';
            subjects = [1002 1005 1006 1007];
        elseif SRType == 2
            levels = [0:.1:1];
            PlotTitle = ' VSR Auditory Task'; xTitle = 'VSR Level, mA'; yTitle = 'Auditory Threshold, dB';
            subjects = [1002 1003 1008];
        else
        end
    case 3 % define information for tactile task
        if SRType == 1
            levels = [0 30:5:80];
            PlotTitle = ' ASR Tactile Task'; xTitle = 'ASR Level, dB'; yTitle = 'Tactile Threshold, Displace Percentage';
            subjects = [1003 1006];
        elseif SRType == 2
            levels = [0:.1:1];
            PlotTitle = ' VSR Tactile Task'; xTitle = 'VSR Level, mA'; yTitle = 'Tactile Threshold, Displace Percentage';
            subjects = [1002 1005];
        else
        end
    case 4 % define information for vestibular task
        if SRType == 1
            levels = [0 30:5:80];
            PlotTitle = ' ASR Vestibular Task'; xTitle = 'ASR Level, dB'; yTitle = 'Vestibular Threshold, deg';
            subjects = [1002 1005];
        elseif SRType == 2
            levels = [0:.1:1];
            PlotTitle = ' VSR Vestibular Task'; xTitle = 'VSR Level, mA'; yTitle = 'Vestibular Threshold, deg';
        else
        end
end

%% Loop through all the subjects and plot their data
currentFolder=pwd;
typeFolder='RA1Data';
for d=1:length(subjects)
    %% define files to be read and define initial threshold guesses
    subID=subjects(d);
    subFold=num2str(subID);
    disp(['Analyzing ' num2str(subID)])
    switch TestType
        case 1 % define visual folder information
            mu_guess = 0.2; sigma_guess = 0.05; N=50;
            if SRType == 1
                parseFold1=['BSRA1_ASR_' num2str(subID) '_Visual_0mA_'];
                parseFold2=['dB.xls'];
            elseif SRType == 2
                parseFold1=['BSRA1_VSR_' num2str(subID) '_Visual_'];
                parseFold2=['mA_0dB.xls'];
            else
            end
        case 2 % define auditory folder information
            mu_guess = 8; sigma_guess = 3; N=100;
            if SRType == 1
                parseFold1=['BSRA1_ASR_' num2str(subID) '_Auditory_0mA_'];
                parseFold2=['dB.xls'];
            elseif SRType == 2
                parseFold1=['BSRA1_VSR_' num2str(subID) '_Auditory_'];
                parseFold2=['mA_0dB.xls'];
            else
            end
        case 3 % define tactile folder information
            mu_guess = 0.3; sigma_guess = 0.1; N=50;
            if SRType == 1
                parseFold1=['BSRA1_ASR_' num2str(subID) '_Tactile_0mA_'];
                parseFold2=['dB.xls'];
            elseif SRType == 2
                parseFold1=['BSRA1_VSR_' num2str(subID) '_Tactile_'];
                parseFold2=['mA_0dB.xls'];
            else
            end
        case 4
            mu_guess = 0; sigma_guess = 0; N=100;
            if SRType == 1
                parseFold1=['BSRA1_ASR_' num2str(subID) '_Vestibular_0mA_'];
                parseFold2=['dB.xls'];
            elseif SRType == 2
                parseFold1=['BSRA1_VSR_' num2str(subID) '_Vestibular_'];
                parseFold2=['mA_0dB.xls'];
            else
            end
    end
    %% Loop through levels, calculate thresholds, simulate oh my!
    guess_rate = 0.5;  % the level at which subject guesses correctly, normally 0.5
    for i=1:length(levels)
        
        % Determine the file you are entering... slightly ugly here
        level=levels(i);
        if sign(level) == -1 % folder definitions for negative SR levels
            fileID=[parseFold1 'neg' num2str(abs(level)) parseFold2];
            fid=fullfile(currentFolder,typeFolder,subFold,fileID);
            master=importdata(fid);
        elseif floor(level) ~= level % definitions for non-integer SR levels
            fileID=[parseFold1 num2str(floor(level)) 'p' num2str(10*rem(level,1)) parseFold2];
            if TestType == 2
                fileID=[parseFold1 num2str(floor(level)) 'p' num2str(10*rem(level,1)) parseFold2];
            end
            fid=fullfile(currentFolder,typeFolder,subFold,fileID);
            master=importdata(fid);
        else % definition for all whole number SR levels
            fileID=[parseFold1 num2str(abs(level)) parseFold2];
            fid=fullfile(currentFolder,typeFolder,subFold,fileID);
            master=importdata(fid);
        end
        
        % separate out the data, do NOT take absolute value for auditory
        % this was already accounted for in the test code for auditory
        if TestType == 2
            X = master.data.Sheet1(:,2);
            cor = master.data.Sheet1(:,4);
        else
            X = master.data.Sheet1(:,2);
            X = abs(X);
            cor = master.data.Sheet1(:,4);
        end
        
        % Fit binary thresholds and find their jack knife bars
        disp(['Data binary, noise is ' num2str(level)])
        x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mu_guess, sigma_guess]);
        mu_bi = x(1);
        sigma_bi = x(2);
        %     thresh_794_bi = icdf('norm', (0.794 - guess_rate)/(1-guess_rate), mu_bi, sigma_bi);
        disp(['Jackknife binary, noise is ' num2str(level)])
        [musi_bi] = jackknife(@just794_bi,X,cor,mu_guess,sigma_guess);
        mu_se_bi = sqrt( (N-1)/N*sum( (mu_bi - musi_bi(:,1)).^2) );
        sigma_se_bi = sqrt( (N-1)/N*sum( (sigma_bi - musi_bi(:,2)).^2 ) );
        % shove and save our data in cell arrays
        AllData{i,1} = X; % stim levels
        AllData{i,2} = cor; % correct?
        AllData{i,3} = [mu_bi sigma_bi]; % mus
        AllData{i,4} = [mu_se_bi sigma_se_bi]; % binary jack knife
    end
    
    %% Parse out AllData for plotting
    
    level_vect=levels;
    for i=1:length(levels)
        % Pool out binary information
        BiThresh=AllData{i,3}; BiJack=AllData{i,4};
        mubi_vect(i)=BiThresh(1); sigmabi_vect(i)=BiThresh(2);
        mubi_jack_vect(i)=BiJack(1); sigmabi_jack_vect(i)=BiJack(2);;
        BiThresh=[]; BiJack=[];
        
    end
    
    %% Develop Vectors and reorganize for plotting, yes it's ugly, but there is a method to the madness
    
    % VECTORS ARE ONLY USED FOR PLOTTING PURPOSES!!!
    if TestType == 2 % auditory reorganizing, thank god we only need to do it for this one
        if SRType == 1
            level_vect(2)=levels(1)-10; % no noise condition, ASR
            level_vect([1 2])=level_vect([2 1]); % we need to swap due to negative condition
            mubi_vect([1 2])=mubi_vect([2 1]); sigmabi_vect([1 2])=sigmabi_vect([2 1]);
            mubi_jack_vect([1 2])=mubi_jack_vect([2 1]); sigmabi_jack_vect([1 2])=sigmabi_jack_vect([2 1]);
        elseif SRType == 2
            level_vect(1)=levels(1)-0.1; % no noise condition, VSR
        end
    else % Other testing reorganizing
        if SRType == 1
            level_vect(1)=levels(1)+20; % no noise condition, ASR
        elseif SRType == 2
            level_vect(1)=levels(1)-0.1; % no noise condition, VSR
        end
    end
    % vectors to become reference dash line
    shamline_bi=mubi_vect(1)*ones(1,length(level_vect(2:end)));
    
    %% Define Boundaries and tick marks
    switch TestType
        case 1
            textscale = -.1;
            if SRType == 1
                tickMarks=[20 30:5:80]; xlimit=[17.5 82.5]; ylimit=[-.15 .6]; Toffset=1;
                xticklabs={'No Noise' '30' '35' '40' '45' '50' '55' '60' '65' '70' '75' '80'};
            elseif SRType == 2
                tickMarks=[-.1 .1:.1:1]; xlimit=[-.3 1.2]; ylimit=[-.15 .6]; Toffset=.025;
                xticklabs={'No Noise' '.1' '.2' '.3' '.4' '.5' '.6' '.7' '.8' '.9' '1'};
            end
        case 2
            textscale = -8;
            if SRType == 1
                tickMarks=[-15 -5:5:40]; xlimit=[-20 45]; ylimit=[-10 20]; Toffset=1;
                xticklabs={'No Noise' '-5' '0' '5' '10' '15' '20' '25' '30' '35' '40'};
            elseif SRType == 2
                tickMarks=[-.1 .1:.1:1]; xlimit=[-.3 1.2]; ylimit=[-10 20]; Toffset=.025;
                xticklabs={'No Noise' '.1' '.2' '.3' '.4' '.5' '.6' '.7' '.8' '.9' '1'};
            end
        case 3
            textscale = -.1;
            if SRType == 1
                tickMarks=[20 30:5:80]; xlimit=[17.5 82.5]; ylimit=[-.15 1]; Toffset=1;
                xticklabs={'No Noise' '30' '35' '40' '45' '50' '55' '60' '65' '70' '75' '80'};
            elseif SRType == 2
                tickMarks=[-.1 .1:.1:1]; xlimit=[-.3 1.2]; ylimit=[-.15 1]; Toffset=.025;
                xticklabs={'No Noise' '.1' '.2' '.3' '.4' '.5' '.6' '.7' '.8' '.9' '1'};
            end
        case 4
            textscale = -.1;
            if SRType == 1
                tickMarks=[20 30:5:80]; xlimit=[17.5 82.5]; ylimit=[-.15 .6]; Toffset=1;
                xticklabs={'No Noise' '30' '35' '40' '45' '50' '55' '60' '65' '70' '75' '80'};
            elseif SRType == 2
                tickMarks=[-.1 .1:.1:1]; xlimit=[-.3 1.2]; ylimit=[-.15 .6]; Toffset=.025;
                xticklabs={'No Noise' '.1' '.2' '.3' '.4' '.5' '.6' '.7' '.8' '.9' '1'};
            end
    end
    
    %% Plot all SR curves
    subplot(2,1,d)
    plot(level_vect,mubi_vect,['*k'],'LineWidth',1.1)
    hold on
    plot(level_vect(2:end),mubi_vect(2:end),['--k'])
    plot(level_vect(2:end),shamline_bi,['--k'],'LineWidth',1.05)
    errorbar(level_vect,mubi_vect,mubi_jack_vect*1.96,'k');
    text(level_vect'-Toffset, textscale*ones(length(level_vect),1),num2str(round(sigmabi_vect',2)))
    %legVect(d,:)=[num2str(subID) ' ' num2str(round(mubi_vect(1),3))];
    hold off
    xticks(tickMarks); xlim(xlimit); ylim(ylimit); xticklabels(xticklabs);
    xlabel(xTitle); ylabel(yTitle); title([num2str(subjects(d)) ' ' PlotTitle ' Binary, Jack'])
    %legend(h1,legVect,'Location','northeastoutside','FontSize',11)
    set(gca,'FontSize',11)
end