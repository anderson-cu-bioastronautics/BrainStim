function plotPreferences2(NoisePrefer_mat,testData_mat,mdlD,numSubs,titCond,saveFolder) % Do the outlier analysis and plotting later

%% Define parameters
yvalues = 1:(numSubs);
titles = {'DSST', 'LOT', 'MPT', 'MRT', 'NBACK', 'PVT', 'VOLT'};
xvalues = {'nGVS', 'AWN', 'MMSR'};
markers = ['o' '+' '*' '.' 'x' 'p' 'h' '^' '>' '<' 'v' 's' 'd'];
ColorVect = ['r','b','m'];
xLims = [-5 5];
%xLims = [0 6];

%% Collapse, plot, and regress the data

    figure()
    plotTitle = ['Collapsed % Difference against ' titCond];
    sgtitle(plotTitle);
    
    % plot feedback information
    % subplot(1,3,1);
    for kay = 1:3
        NoisePrefer = NoisePrefer_mat{kay}; testData = testData_mat{kay}; mdl = mdlD{kay};
        y1(kay) = plot(NoisePrefer,testData,'*','color',ColorVect(kay),'Marker', markers(kay))
        hold on
        xFit = linspace(min(NoisePrefer), max(NoisePrefer), 1000);
        % Get the estimated yFit value for each of those 1000 new x locations.
        yFit = polyval([mdl.Coefficients{2,1} mdl.Coefficients{1,1}],xFit);
        plot(xFit,yFit,'-','color',ColorVect(kay),'LineWidth',1.35); % Plot fitted line.
        NoisePrefer = []; FeedData = []; mdl = []; xFit = []; yFit = [];
    end
    hold off
    legend(y1,'nGVS','AWN,','MMSR')
    % title('Feedback');
    ylabel('Normalized Percentage Improvement'); xlabel('Preference Score');
    xlim(xLims)
    set(gca,'FontSize',13)
    
    figName = ['AllPreference' [titCond] '.jpg'];
    set(gcf, 'PaperUnits', 'inches');
   % y_width=7 ;x_width=10;
   % set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
   saveas(gcf, fullfile(saveFolder, figName));
    