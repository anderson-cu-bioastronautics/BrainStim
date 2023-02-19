% function to provide plot data of all our information

function PlotData(AllData,subID,TestType,SRType,levels,PlotTitle,xTitle,yTitle,plotColor)
colors;
%% Parse out AllData

level_vect=levels;
for i=1:length(levels)
    % Pool out binary information
    BiThresh=AllData{i,3}; BiJack=AllData{i,5}; BiBoot=AllData{i,7};
    mubi_vect(i)=BiThresh(1); sigmabi_vect(i)=BiThresh(2);
    mubi_jack_vect(i)=BiJack(1); sigmabi_jack_vect(i)=BiJack(2);
    mubi_boot_vect(i,1)=BiBoot(1); mubi_boot_vect(i,2)=BiBoot(2); sigmabi_boot_vect(i,1)=BiBoot(3); sigmabi_boot_vect(i,2)=BiBoot(4);
    BiThresh=[]; BiJack=[]; BiBoot=[];
    
    % Pool out lapse information
    HatThresh=AllData{i,4}; HatJack=AllData{i,6}; HatBoot=AllData{i,8};
    muhat_vect(i)=HatThresh(1); sigmahat_vect(i)=HatThresh(2);
    muhat_jack_vect(i)=HatJack(1); sigmahat_jack_vect(i)=HatJack(2);
    muhat_boot_vect(i,1)=HatBoot(1); muhat_boot_vect(i,2)=HatBoot(2); sigmahat_boot_vect(i,1)=HatBoot(3); sigmahat_boot_vect(i,2)=HatBoot(4);
    HatThresh=[]; HatJack=[]; HatBoot=[];
end

%% Develop Vectors and reorganize for plotting, yes it's ugly, but there is a method to the madness

% VECTORS ARE ONLY USED FOR PLOTTING PURPOSES!!!
if TestType == 2 % auditory reorganizing, thank god we only need to do it for this one
    if SRType == 1
        level_vect(2)=levels(1)-10; % no noise condition, ASR
        level_vect([1 2])=level_vect([2 1]); % we need to swap due to negative condition
        mubi_vect([1 2])=mubi_vect([2 1]); sigmabi_vect([1 2])=sigmabi_vect([2 1]);
        mubi_jack_vect([1 2])=mubi_jack_vect([2 1]); sigmabi_jack_vect([1 2])=sigmabi_jack_vect([2 1]);
        mubi_boot_vect([1 2],:)=mubi_boot_vect([2 1],:); sigmabi_boot_vect([1 2],:)=sigmabi_boot_vect([2 1],:);
        muhat_vect([1 2])=muhat_vect([2 1]); sigmahat_vect([1 2])=sigmahat_vect([2 1]);
        muhat_jack_vect([1 2])=muhat_jack_vect([2 1]); sigmahat_jack_vect([1 2])=sigmahat_jack_vect([2 1]);
        muhat_boot_vect([1 2],:)=muhat_boot_vect([2 1],:); sigmahat_boot_vect([1 2],:)=sigmahat_boot_vect([2 1],:);
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
shamline_hat=muhat_vect(1)*ones(1,length(level_vect(2:end)));

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

%% Plot the SR curves (in subplots)
figure
hold on;
% Plot the confidence interval of the baseline value across the entire
% figure
rect_x_start = level_vect(1);
rect_y_start = mubi_vect(1) - mubi_jack_vect(1)*1.96;
rect_width = level_vect(length(level_vect)) - rect_x_start;
rect_height = mubi_jack_vect(1)*3.92;
rect_pos = [rect_x_start,rect_y_start,rect_width,rect_height];
% Draw the rectangle
rectangle('Position',rect_pos,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]);

if plotColor ==0
    %subplot(2,2,1) % binary fits with jack knife error bars
    plot(level_vect,mubi_vect,'*k','LineWidth',1.1)
    plot(level_vect(2:end),mubi_vect(2:end),'--','Color',[.17 .17 .17])
    plot(level_vect(2:end),shamline_bi,'--k','LineWidth',1.05)
    errorbar(level_vect,mubi_vect,mubi_jack_vect*1.96,'k')
    text(level_vect'-Toffset, textscale*ones(length(level_vect),1),num2str(round(sigmabi_vect',2)))
    xticks(tickMarks); xlim(xlimit); ylim(ylimit); xticklabels(xticklabs);
    xlabel(xTitle); ylabel(yTitle); title([num2str(subID) PlotTitle ' Binary, Jack'])
        set(gca,'FontSize',12.5)
    %     set(gca,'color',[250,245,252]/255)
    %     set(gcf,'Color',[250,245,252]/255)
    hold off
else
    %subplot(2,2,1) % binary fits with jack knife error bars
    plot(level_vect,mubi_vect,'*','Color',colorOrder(plotColor,:)/255,'LineWidth',1.1)
    plot(level_vect(2:end),mubi_vect(2:end),'--','Color',[.17 .17 .17])
    plot(level_vect(2:end),shamline_bi,'--','Color',colorOrder(plotColor,:)/255,'LineWidth',1.05)
    errorbar(level_vect,mubi_vect,mubi_jack_vect*1.96,'Color',colorOrder(plotColor,:)/255)
    text(level_vect'-Toffset, textscale*ones(length(level_vect),1),num2str(round(sigmabi_vect',2)))
    xticks(tickMarks); xlim(xlimit); ylim(ylimit); xticklabels(xticklabs);
    xlabel(xTitle); ylabel(yTitle); title([num2str(subID) PlotTitle ' Binary, Jack'])
    set(gca,'FontSize',12.5)
    %     set(gca,'color',[250,245,252]/255)
    %     set(gcf,'Color',[250,245,252]/255)
    hold off
end

[lowest_val,lowest_ind]=min(mubi_vect(2:end));
lowest_level=level_vect(lowest_ind+1);
disp(['The retest noise level is ' num2str(lowest_level)]);
% subplot(2,2,2) % lapse fits with jack knife error bars
% plot(level_vect,muhat_vect,'*k','LineWidth',1.1)
% hold on
% plot(level_vect(2:end),muhat_vect(2:end),'--','Color',[.17 .17 .17])
% plot(level_vect(2:end),shamline_hat,'--k','LineWidth',1.05)
% errorbar(level_vect,muhat_vect,muhat_jack_vect*1.96,'k')
% text(level_vect'-Toffset, textscale*ones(length(level_vect),1),num2str(round(sigmahat_vect',2)))
% xticks(tickMarks); xlim(xlimit); ylim(ylimit); xticklabels(xticklabs);
% xlabel(xTitle); ylabel(yTitle); title([num2str(subID) PlotTitle ' Lapse, Jack'])
% hold off
%
% subplot(2,2,3) % binary fits with bootstrap error bars
% plot(level_vect,mubi_vect,'*k','LineWidth',1.1)
% hold on
% plot(level_vect(2:end),mubi_vect(2:end),'--','Color',[.17 .17 .17])
% plot(level_vect(2:end),shamline_bi,'--k','LineWidth',1.05)
% errorbar(level_vect,mubi_vect,mubi_vect-mubi_boot_vect(:,1)',mubi_boot_vect(:,2)'-mubi_vect,'k')
% text(level_vect'-Toffset, textscale*ones(length(level_vect),1),num2str(round(sigmabi_vect',2)))
% xticks(tickMarks); xlim(xlimit); ylim(ylimit); xticklabels(xticklabs);
% xlabel(xTitle); ylabel(yTitle); title([num2str(subID) PlotTitle ' Binary, Boot'])
% hold off
%
% subplot(2,2,4) % lapse fits with bootstrap error bars
% plot(level_vect,muhat_vect,'*k','LineWidth',1.1)
% hold on
% plot(level_vect(2:end),muhat_vect(2:end),'--','Color',[.17 .17 .17])
% plot(level_vect(2:end),shamline_hat,'--k','LineWidth',1.05)
% errorbar(level_vect,muhat_vect,muhat_vect-muhat_boot_vect(:,1)',muhat_boot_vect(:,2)'-muhat_vect,'k')
% text(level_vect'-Toffset, textscale*ones(length(level_vect),1),num2str(round(sigmahat_vect',2)))
% xticks(tickMarks); xlim(xlimit); ylim(ylimit); xticklabels(xticklabs);
% xlabel(xTitle); ylabel(yTitle); title([num2str(subID) PlotTitle ' Lapse, Boot'])
% hold off
% %set(gca,'FontSize',12)

%% Plot the psychometric curves
figure


% Sequence intervals
seqInt = floor(100/length(levels));

% Determine the plot color sequence
switch plotColor
    % Black
    case 0
        seqOrderInit = transpose(linspace(1, 0, length(levels)));
        seqOrder = [seqOrderInit,seqOrderInit,seqOrderInit];
    % Red
    case 1
        seqOrder = redseq(1:seqInt:end);
    % Blue
    case 2
        seqOrder = blueseq(1:seqInt:end);
    % Green
    case 3
        seqOrder = greenseq(1:seqInt:end);
    % Purple
    case 4
        seqOrder = purpleseq(1:seqInt:end);
    % Lavender
    case 5
        seqOrder = lavseq(1:seqInt:end);
end
subplot(1,2,1)
for i=1:length(levels)
    plot(AllData{i,9},AllData{i,10},'Color',seqOrder(i,:));
    hold on
end
xlabel('Stimulus'); ylabel('Likelihood of Correct Response');
%xlim(xlimit);
ylim([.45 1.005]); title([num2str(subID) PlotTitle ' Binary Curves'])

subplot(1,2,2)
for i=1:length(levels)
    plot(AllData{i,9},AllData{i,11},'Color',seqOrder(i,:))
    hold on
end
xlabel('Stimulus'); ylabel('Likelihood of Correct Response');
%xlim(xlimit);
ylim([.45 1.005]); title([num2str(subID) PlotTitle ' Lapse Curves'])
%set(gca,'FontSize',12)