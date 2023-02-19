%% SR Curve Fits Take 2
%  I swear to god if this computer destroys another flash drive I'll have a
%  conniption

%  Written by Rachel Rise
%  Orignially written sometime in March but it's April 16, 2020 now
%  Lmao now it's May 26

%% Let's go

tic;
clc;
clear;
close all;
format longg;

options = optimoptions('fmincon','Display','off');
warning('off');
% Matrices for all the fitted values
ASR_vis = [];
ASR_vis_guesses = [];
ASR_vis_closeParams = [];
ASR_aud = [];
ASR_aud_guesses = [];
ASR_aud_closeParams = [];
ASR_tac = [];
ASR_tac_guesses = [];
ASR_tac_closeParams = [];
ASR_vst = [];
ASR_vst_guesses = [];
ASR_vst_closeParams = [];
VSR_vis = [];
VSR_vis_guesses = [];
VSR_vis_closeParams = [];
VSR_aud = [];
VSR_aud_guesses = [];
VSR_aud_closeParams = [];
VSR_tac = [];
VSR_tac_guesses = [];
VSR_tac_closeParams = [];
VSR_vst = [];
VSR_vst_guesses = [];
VSR_vst_closeParams = [];

% Identify which variables are global
global levels thresholds TestType SRType method widthThreshold1 widthThreshold2

% What method of determining curve width are we using?
% Options are:
%    None: no constraints on curve width or depth
%    Derivative: uses the minimum and maximum slope of the SR curve to
%    estimate the width
%    Percentage: uses the points where the curve crosses certain heights in
%    relation to dip depth
%    Integral: integrates the area under the curve and divides it by the
%    dip depth to determine the width.
method = 'Derivative';

% If using the 'Percentile' method, choose at what percentage of the dip
% depth to define curve width. Threshold 1 is the point where the curve
% dips down, and threshold 2 is where it comes back up.
widthThreshold1 = 0.9;
widthThreshold2 = 0.9;

% This is definitely not efficient. Here's a matrix of all the things we
% want to store and analyze:
%    [1]  Subject #
%    [2]  SR Type
%    [3]  Test Type
%    [4]  A0 guess
%    [5]  lambda guess
%    [6]  w0 guess
%    [7]  m guess
%    [8]  s guess
%    [9]  f guess
%    [10] B guess
%    [11] A0 fit
%    [12] lambda fit
%    [13] w0 fit
%    [14] m fit
%    [15] s fit
%    [16] f fit
%    [17] B fit
%    [18] Cost J
%    [19] Exit condition flag
%             All algorithms:
%               1  First order optimality conditions satisfied.
%               0  Too many function evaluations or iterations.
%              -1  Stopped by output/plot function.
%              -2  No feasible point found.
%             Trust-region-reflective, interior-point, and sqp:
%               2  Change in X too small.
%             Trust-region-reflective:
%               3  Change in objective function too small.
%             Active-set only:
%               4  Computed search direction too small.
%               5  Predicted change in objective function too small.
%             Interior-point and sqp:
%              -3  Problem seems unbounded.
%    [20] Closeness to minimum J flag
%               3  Exactly equal to minimum value
%               2  Within 99% of minimum value
%               1  Within 95% of minimum value
%               0  Outside 95% of minimum
%    [21] Closeness to best parameters
%               0  Exactly equal to best parameters
%               1  Not equal to best parameters

% Horribly inefficient, but I can't think of a better way to organize all
% this
fitTrackingMatrix = [];


%% Let's loop through the SR Types

% Let's set up looooooooops

for SRType = 1:2
    
    % And the test types
    for TestType = 1:4
        
        % And the subjects!
        for subID = 1:13
            
            % Determine what SR type we're going to label our plots with
            switch(SRType)
                
                case 1 % Auditory
                    
                    % What type of SR file are we looking for?
                    SR_string = 'ASR';
                    
                case 2 % Vestibular
                    
                    % What type of SR file are we looking for?
                    SR_string = 'VSR';
                    
            end
            
            % Alright, let's load the subject data.
            loadDataStr = [num2str(subID),'_',SR_string,'_',num2str(TestType)];
            
            % Don't load it if the file doesn't exist
            if exist([loadDataStr,'.mat'], 'file') == 0
                
                % File does not exist, skip this iteration of the loop
                continue;
            else
                
                % The file exists! Load it!
                load(loadDataStr);
                
                % Extract the data we need
                dataSize = size(AllData);
                n = dataSize(1);
                levels = zeros(1,n);
                thresholds = zeros(1,n);
                
                % Populate the data
                for i = 1:n
                    
                    data = AllData{i,3};
                    
                    levels(i) = data(3);
                    
                    thresholds(i) = data(1);
                    
                end
                
                
            end
            
            % Provide some guesses for our curve fit
            switch SRType
                
                case 1
                    
                    switch TestType
                        
                        case 1
                            
                            % The guesses (ASR Visual)
                            A0_guesses = [1 10 100 1000];
                            lambda_guesses = [1 5 10 50];
                            w0_guesses = [1 5 10 50];
                            f_guesses = [30.01 40.01 50.01 60.01];
                            B_guess = thresholds(1);
                            
                            % The limits on the cost function
                            LowerBounds = [1e-6;
                                1;
                                1;
                                30;
                                0];
                            
                            UpperBounds = [Inf;
                                Inf;
                                Inf;
                                80;
                                1];
                            
                            
                            
                            % Title, labels, and axes
                            titlestr = ['Subject ',num2str(subID),' ASR Visual'];
                            xstr = 'ASR Level (dB)';
                            ystr = 'Threshold Contrast Level';
                            axes = [-5 85 0 1];
                            
                        case 2
                            % The guesses (ASR Auditory)
                            A0_guesses = [1 10 100 1000];
                            lambda_guesses = [1 5 10 50];
                            w0_guesses = [1 5 10 50];
                            m_guess = 0;
                            s_guess = 30;
                            f_guesses = [1e-6 5.1 10.1 15.1 20.1 25.1 30.1];
                            B_guess = 8;
                            
                            % The limits on the cost function
                            LowerBounds = [1e-6;
                                1;
                                1;
                                0;
                                10;
                                -5;
                                -10];
                            
                            UpperBounds = [Inf;
                                Inf;
                                Inf;
                                1;
                                40;
                                40;
                                25];
                            
                            % Change the first level so the curve fit works
                            % right
                            levels(1) = -15;
                            
                            % Title, labels, and axes
                            titlestr = ['Subject ',num2str(subID),' ASR Auditory'];
                            xstr = 'ASR Level (dB)';
                            ystr = 'Threshold (dB)';
                            axes = [-20 45 -10 40];
                            
                        case 3
                            % The guesses (ASR Tactile)
                            A0_guesses = [1 10 100 1000];
                            lambda_guesses = [1 5 10 50];
                            w0_guesses = [1 5 10 50];
                            f_guesses = [30.1 40.1 50.1 60.1];
                            B_guess = 0.2;
                            
                            % The limits on the cost function
                            LowerBounds = [1e-6;
                                1;
                                1;
                                30;
                                0];
                            
                            UpperBounds = [Inf;
                                Inf;
                                Inf;
                                80;
                                1];
                            
                            % Title, labels, and axes
                            titlestr = ['Subject ',num2str(subID),' ASR Tactile'];
                            xstr = 'ASR Level (dB)';
                            ystr = 'Threshold Normalized Amplitude';
                            axes = [-5 85 0 1];
                            
                        case 4
                            % The guesses (ASR Vestibular)
                            A0_guesses = 1;
                            lambda_guesses = 5;
                            w0_guesses = 5;
                            f_guesses = 10;
                            B_guess = 8;
                            
                            % The limits on the cost function
                            LowerBounds = [1e-6;
                                1;
                                1;
                                30;
                                0];
                            
                            UpperBounds = [Inf;
                                Inf;
                                Inf;
                                80;
                                1];
                            
                            % Title, labels, and axes
                            titlestr = ['Subject ',num2str(subID),' ASR Vestibular'];
                            xstr = 'ASR Level (dB)';
                            ystr = 'Threshold Sensitivity (deg/s)';
                            axes = [-10 10 0 10];
                    end
                    
                case 2
                    
                    switch TestType
                        
                        case 1
                            % The guesses (VSR Visual)
                            A0_guesses = [1 10 100 1000];
                            lambda_guesses = [1 5 10 50];
                            w0_guesses = [1 5 10 50];
                            f_guesses = [0.01 0.11 0.21 0.31 0.41 0.51 0.61 0.71 0.81];
                            B_guess = 0.2;
                            
                            % The limits on the cost function
                            LowerBounds = [1e-6;
                                1;
                                1;
                                0;
                                0];
                            
                            UpperBounds = [Inf;
                                Inf;
                                Inf;
                                1;
                                1];
                            
                            
                            % Title, labels, and axes
                            titlestr = ['Subject ',num2str(subID),' VSR Visual'];
                            xstr = 'VSR Level (mA)';
                            ystr = 'Threshold Contrast Level';
                            axes = [-0.1 1.1 0 1];
                            
                        case 2
                            % The guesses (VSR Auditory)
                            A0_guesses = [1 10 100 1000];
                            lambda_guesses = [1 5 10 50];
                            w0_guesses = [1 5 10 50];
                            f_guesses = [0.01 0.11 0.21 0.31 0.41 0.51 0.61 0.71 0.81];
                            B_guess = 8;
                            
                            % The limits on the cost function
                            LowerBounds = [1e-6;
                                1;
                                1;
                                0;
                                0];
                            
                            UpperBounds = [Inf;
                                Inf;
                                Inf;
                                1;
                                40];
                            
                            
                            
                            % Title, labels, and axes
                            titlestr = ['Subject ',num2str(subID),' VSR Auditory'];
                            xstr = 'VSR Level (mA)';
                            ystr = 'Threshold (dB)';
                            axes = [-0.1 1.1 0 40];
                            
                        case 3
                            % The guesses (VSR Tactile)
                            A0_guesses = 1;
                            lambda_guesses = 5;
                            w0_guesses = 5;
                            f_guesses = 0.3;
                            B_guess = 0.2;
                            
                            % The limits on the cost function
                            LowerBounds = [1e-6;
                                1;
                                1;
                                0;
                                0];
                            
                            UpperBounds = [Inf;
                                Inf;
                                Inf;
                                1;
                                1];
                            
                            % Title, labels, and axes
                            titlestr = ['Subject ',num2str(subID),' VSR Visual'];
                            xstr = 'VSR Level (mA)';
                            ystr = 'Threshold Normalized Amplitude';
                            axes = [-0.1 1.1 0 1];
                            
                        case 4
                            % The guesses (VSR Vestibular)
                            A0_guesses = [1 10 100 1000];
                            lambda_guesses = [1 5 10 50];
                            w0_guesses = [1 5 10 50];
                            f_guesses = [0 0.2 0.4 0.6];
                            B_guess = 8;
                            
                            % The limits on the cost function
                            LowerBounds = [1e-6;
                                1;
                                1;
                                30;
                                0];
                            
                            UpperBounds = [Inf;
                                Inf;
                                Inf;
                                80;
                                1];
                            
                            % Title, labels, and axes
                            titlestr = ['Subject ',num2str(subID),' VSR Vestibular'];
                            xstr = 'VSR Level (mA)';
                            ystr = 'Threshold Sensitivity (deg/s)';
                            axes = [-0.1 1.1 0 1];
                    end
            end
            
            % How many fit iterations are we going to do?
            fitIterations = length(A0_guesses)*length(lambda_guesses)*length(w0_guesses)*length(f_guesses);
            
            % Allocate a fit tracking matrix in this loop
            fitTrackingMatrixLocal = zeros(20,fitIterations);
            
            % Set up indexing for each fit iteration
            fitIndex = 1;
            
            % Loop through the guesses and fit the curve with each one
            for i = 1:length(A0_guesses)
                
                A0_guess = A0_guesses(i);
                
                for j = 1:length(lambda_guesses)
                    
                    lambda_guess = lambda_guesses(j);
                    
                    for k = 1:length(w0_guesses)
                        
                        w0_guess = w0_guesses(k);
                        
                        for l = 1:length(f_guesses)
                            
                            f_guess = f_guesses(l);  
                            
                            % Fit the curve!
                            
                            % For ASR auditory, we include masking
                            if (TestType == 2 && SRType == 1)
                                % Compile the guesses
                                guesses = [A0_guess lambda_guess w0_guess m_guess s_guess f_guess B_guess];
                                [params, J, flag] = fmincon(@(X) SR_curve_fit_freeB(X, levels, thresholds), guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth,options);
                                
                                % Allocate everything to the local fit
                                % tracking matrix
                                fitTrackingMatrixLocal(1,fitIndex) = subID;
                                fitTrackingMatrixLocal(2,fitIndex) = SRType;
                                fitTrackingMatrixLocal(3,fitIndex) = TestType;
                                fitTrackingMatrixLocal(4,fitIndex) = A0_guess;
                                fitTrackingMatrixLocal(5,fitIndex) = lambda_guess;
                                fitTrackingMatrixLocal(6,fitIndex) = w0_guess;
                                fitTrackingMatrixLocal(7,fitIndex) = m_guess;
                                fitTrackingMatrixLocal(8,fitIndex) = s_guess;
                                fitTrackingMatrixLocal(9,fitIndex) = f_guess;
                                fitTrackingMatrixLocal(10,fitIndex) = B_guess;
                                fitTrackingMatrixLocal(11,fitIndex) = params(1);
                                fitTrackingMatrixLocal(12,fitIndex) = params(2);
                                fitTrackingMatrixLocal(13,fitIndex) = params(3);
                                fitTrackingMatrixLocal(14,fitIndex) = params(4);
                                fitTrackingMatrixLocal(15,fitIndex) = params(5);
                                fitTrackingMatrixLocal(16,fitIndex) = params(6);
                                fitTrackingMatrixLocal(17,fitIndex) = params(7);
                                fitTrackingMatrixLocal(18,fitIndex) = J;
                                fitTrackingMatrixLocal(19,fitIndex) = flag;
                                
                                % For other conditions, don't include masking
                            else
                                
                                guesses = [A0_guess lambda_guess w0_guess f_guess B_guess];
                                [params, J, flag] = fmincon(@(X) SR_curve_fit_nomask(X, levels, thresholds), guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth,options);
                                
                                % Allocate everything to the local fit
                                % tracking matrix
                                fitTrackingMatrixLocal(1,fitIndex) = subID;
                                fitTrackingMatrixLocal(2,fitIndex) = SRType;
                                fitTrackingMatrixLocal(3,fitIndex) = TestType;
                                fitTrackingMatrixLocal(4,fitIndex) = A0_guess;
                                fitTrackingMatrixLocal(5,fitIndex) = lambda_guess;
                                fitTrackingMatrixLocal(6,fitIndex) = w0_guess;
                                fitTrackingMatrixLocal(7,fitIndex) = 0;
                                fitTrackingMatrixLocal(8,fitIndex) = 0;
                                fitTrackingMatrixLocal(9,fitIndex) = f_guess;
                                fitTrackingMatrixLocal(10,fitIndex) = B_guess;
                                fitTrackingMatrixLocal(11,fitIndex) = params(1);
                                fitTrackingMatrixLocal(12,fitIndex) = params(2);
                                fitTrackingMatrixLocal(13,fitIndex) = params(3);
                                fitTrackingMatrixLocal(14,fitIndex) = 0;
                                fitTrackingMatrixLocal(15,fitIndex) = 0;
                                fitTrackingMatrixLocal(16,fitIndex) = params(4);
                                fitTrackingMatrixLocal(17,fitIndex) = params(5);
                                fitTrackingMatrixLocal(18,fitIndex) = J;
                                fitTrackingMatrixLocal(19,fitIndex) = flag;
                                
                            end
                            
                            % Increment the fit index
                            fitIndex = fitIndex + 1;
                            
                        end
                    end
                end
            end
            
            
            % The minimum cost achieved in by all those guesses
            J_vector = fitTrackingMatrixLocal(18,:);
            Jmin = min(J_vector);
            
            % Find which sets of guesses produced cost functions close to
            % the minimum
            J_closeness = abs(J_vector - Jmin)./Jmin;
            
            % Parameter closeness flag
            paramsMatrix = fitTrackingMatrixLocal(11:17,:);
            bestParams = paramsMatrix(:,min(find(J_vector == Jmin)));
            
            % Evaluate which parameters are actually best
            A0 = bestParams(1);
            lambda = bestParams(2);
            w0 = bestParams(3);
            m = bestParams(4);
            s = bestParams(5);
            f = bestParams(6);
            B = bestParams(7);
            
            % Find which guesses resulted in solutions very close to the
            % best
            guess_closeness_vector = zeros(1,fitIterations);
            
            all_A0 = paramsMatrix(1,:);
            A0_closeness = abs(all_A0 - A0)./A0;
            A0_closeness_flag = (A0_closeness <= 0.01);
            
            all_lambda = paramsMatrix(2,:);
            lambda_closeness = abs(all_lambda - lambda)./lambda;
            lambda_closeness_flag = (lambda_closeness <= 0.01);
            

            all_w0 = paramsMatrix(3,:);
            w0_closeness = abs(all_w0 - w0)./w0;
            w0_closeness_flag = (w0_closeness <= 0.01);
            
            if TestType == 2 && SRType == 1
            all_m = paramsMatrix(4,:);
            m_closeness = abs(all_m - m)./m;
            m_closeness_flag = (m_closeness <= 0.01);
            
            all_s = paramsMatrix(5,:);
            s_closeness = abs(all_s - s)./s;
            s_closeness_flag = (s_closeness <= 0.01);
            else
                m_closeness_flag = ones(1,fitIterations);
                s_closeness_flag = ones(1,fitIterations);
            end
            all_f = paramsMatrix(6,:);
            f_closeness = abs(all_f - f)./f;
            f_closeness_flag = (f_closeness <= 0.01);
            
            all_B = paramsMatrix(7,:);
            B_closeness = abs(all_B - B)./B;
            B_closeness_flag = (B_closeness <= 0.01);
            
            paramsCloseness = A0_closeness_flag.*lambda_closeness_flag.*w0_closeness_flag.*m_closeness_flag.*s_closeness_flag.*f_closeness_flag.*B_closeness_flag;
            
            J_flag = J_closeness <= 0.01;
            
            flag_vector = (fitTrackingMatrixLocal(19,:) ~= 1);
            
            % Fit closeness vector
            fit_closeness_vector = zeros(1,fitIterations);
                
            meep = (~J_flag);
            
            fit_closeness_vector(find(meep)) = 1;
            
            morp = (J_flag).*(~paramsCloseness);
            
            fit_closeness_vector(find(morp)) = 2;
            
            blah = (J_flag).*(paramsCloseness);
            
            fit_closeness_vector(find(blah)) = 3;
            
            fit_closeness_vector(find(J_vector == Jmin)) = 4;
            
            
            fitTrackingMatrixLocal(20,:) = fit_closeness_vector;
            
            
            % the best guesses
            guessMatrix = fitTrackingMatrixLocal(4:10,:);
            bestGuesses = guessMatrix(:,min(find(J_vector == Jmin)));

            if (SRType == 1 && TestType == 2)
                
                disp(titlestr);
                fprintf('   A0:     %16.3f | Guess:   %16.3f\n',A0,bestGuesses(1));
                fprintf('   Lambda: %16.3f | Guess:   %16.3f\n',lambda,bestGuesses(2));
                fprintf('   w0:     %16.3f | Guess:   %16.3f\n',w0,bestGuesses(3));
                fprintf('   m:      %16.3f | Guess:   %16.3f\n',m,bestGuesses(4));
                fprintf('   s:      %16.3f | Guess:   %16.3f\n',s,bestGuesses(5));
                fprintf('   f:      %16.3f | Guess:   %16.3f\n',f,bestGuesses(6));
                fprintf('   B:      %16.3f | Guess:   %16.3f\n',B,bestGuesses(7));
                fprintf('   Cost J: %16.3f\n\n',Jmin);
                fprintf('\n\n');
                
            else
                
                disp(titlestr);
                fprintf('   A0:     %16.3f | Guess:   %16.3f\n',A0,bestGuesses(1));
                fprintf('   Lambda: %16.3f | Guess:   %16.3f\n',lambda,bestGuesses(2));
                fprintf('   w0:     %16.3f | Guess:   %16.3f\n',w0,bestGuesses(3));
                fprintf('   f:      %16.3f | Guess:   %16.3f\n',f,bestGuesses(6));
                fprintf('   B:      %16.3f | Guess:   %16.3f\n',B,bestGuesses(7));
                fprintf('   Cost J: %16.3f\n',Jmin);                
                fprintf('\n\n');
            end
            J_str = sprintf(', J = %4.3f',Jmin);
            
            % Add the fitted values to the appropriate matrix
            switch SRType
                case 1
                    switch TestType
                        case 1
                            ASR_vis = horzcat(ASR_vis, bestParams);
                            ASR_vis_guesses = horzcat(ASR_vis_guesses, bestGuesses);
                            ASR_vis_closeParams = vertcat(ASR_vis_closeParams, fit_closeness_vector);
                        case 2
                            ASR_aud = horzcat(ASR_aud, bestParams);
                            ASR_aud_guesses = horzcat(ASR_aud_guesses, bestGuesses);
                            ASR_aud_closeParams = vertcat(ASR_aud_closeParams, fit_closeness_vector);
                        case 3
                            ASR_tac = horzcat(ASR_tac, bestParams);
                            ASR_tac_guesses = horzcat(ASR_tac_guesses, bestGuesses);
                            ASR_tac_closeParams = vertcat(ASR_tac_closeParams, fit_closeness_vector);
                        case 4
                            ASR_vst = vertcat(ASR_vst, bestParams);
                            ASR_vst_guesses = horzcat(ASR_vst_guesses, bestGuesses);
                            ASR_vst_closeParams = vertcat(ASR_vst_closeParams, fit_closeness_vector);
                    end
                case 2
                    switch TestType
                        case 1
                            VSR_vis = horzcat(VSR_vis, bestParams);
                            VSR_vis_guesses = horzcat(VSR_vis_guesses, bestGuesses);
                            VSR_vis_closeParams = vertcat(VSR_vis_closeParams, fit_closeness_vector);
                            
                        case 2
                            VSR_aud = horzcat(VSR_aud, bestParams);
                            VSR_aud_guesses = horzcat(VSR_aud_guesses, bestGuesses);
                            VSR_aud_closeParams = vertcat(VSR_aud_closeParams, fit_closeness_vector);
                        case 3
                            VSR_tac = horzcat(VSR_tac, bestParams);
                            VSR_tac_guesses = horzcat(VSR_tac_guesses, bestGuesses);
                            VSR_tac_closeParams = vertcat(VSR_tac_closeParams, fit_closeness_vector);
                        case 4
                            VSR_vst = horzcat(VSR_vst, bestParams);
                            VSR_vst_guesses = horzcat(VSR_vst_guesses, bestGuesses);
                            VSR_vst_closeParams = vertcat(VSR_vst_closeParams, fit_closeness_vector);
                    end
            end
            

            % Concatenate the fit tracking matrix
            fitTrackingMatrix = horzcat(fitTrackingMatrix,fitTrackingMatrixLocal);
            
            % Make an x vector
            X = linspace(min(levels), max(levels),1000);
            
            % Generate the fitted curve
            r = (lambda./(sqrt(2)*pi)).*exp(-lambda^2./(2*((X-f).^2)));
            F = B - (A0.*lambda./(X - f).^2).*r./sqrt(4.*r.^2 + w0^2).*(X-f).*(X >=f);
            if(SRType == 1 && TestType == 2)
                mask = m*(X-s).*(X >= s);
                F = F + mask;
            end
            
            % Plot it!
            figure;
            hold on;
            plot(levels,thresholds,'k*','MarkerSize',10,'LineStyle','None');
            plot(X,F,'r','LineWidth',1.5);
            title([titlestr,J_str],'FontSize',20);
            xlabel(xstr,'FontSize',18);
            ylabel(ystr,'FontSize',20);
            axis(axes);
            set(gca,'FontSize',16);
            legend('Data Points','Fitted Curve');
            
            % Create a text box annotation
            if (SRType == 1 && TestType == 2)
                textstr = {sprintf('A0:       %4.3f',A0),
                    sprintf('lambda:   %4.3f',lambda),
                    sprintf('w0:       %4.3f',w0),
                    sprintf('f:        %4.3f',f),
                    sprintf('B:        %4.3f',B),
                    sprintf('m:        %4.3f',m),
                    sprintf('s:        %4.3f',s)};
            else
                textstr = {sprintf('A0:       %4.3f',A0),
                    sprintf('lambda:   %4.3f',lambda),
                    sprintf('w0:       %4.3f',w0),
                    sprintf('f:        %4.3f',f),
                    sprintf('B:        %4.3f',B)};
            end
            
            % Annotation dimensions
            dim = [0.15 0.4 0.3 0.4];
            
            % Make the annotation
            annotation('textbox',dim,'String',textstr,'FontName','FixedWidth','FontSize',14,'FitBoxToText','On','EdgeColor',[1 1 1],'BackgroundColor','none');
            
        end
        
        switch SRType
                case 1
                    switch TestType
                        case 1
                            ASR_vis_guessTracking = vertcat(guessMatrix,ASR_vis_closeParams);
                        case 2
                            ASR_aud_guessTracking = vertcat(guessMatrix,ASR_aud_closeParams);
                        case 3
                            ASR_tac_guessTracking = vertcat(guessMatrix,ASR_tac_closeParams);
                        case 4
                            ASR_vst_guessTracking = vertcat(guessMatrix,ASR_vst_closeParams);
                    end
                case 2
                    switch TestType
                        case 1
                            VSR_vis_guessTracking = vertcat(guessMatrix,VSR_vis_closeParams);
                            
                        case 2
                            VSR_aud_guessTracking = vertcat(guessMatrix,VSR_aud_closeParams);
                        case 3
                            VSR_tac_guessTracking = vertcat(guessMatrix,VSR_tac_closeParams);
                        case 4
                            VSR_vst_guessTracking = vertcat(guessMatrix,VSR_vst_closeParams);
                    end
            end
    end
end
SR = fitTrackingMatrix(2,:);
test = fitTrackingMatrix(3,:);
ASR_vis_indices = find(((SR == 1).*(test == 1)) == 1);
ASR_aud_indices = find(((SR == 1).*(test == 2)) == 1);
ASR_tac_indices = find(((SR == 1).*(test == 3)) == 1);
VSR_vis_indices = find(((SR == 2).*(test == 1)) == 1);
VSR_aud_indices = find(((SR == 2).*(test == 2)) == 1);
VSR_tac_indices = find(((SR == 2).*(test == 3)) == 1);

% ASR visual
figure;
subplot(2,2,1);
loglog(fitTrackingMatrix(4,ASR_vis_indices),fitTrackingMatrix(11,ASR_vis_indices),'.','LineStyle','None','MarkerSize',10);
title('ASR Visual A_0 Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,2);
loglog(fitTrackingMatrix(5,ASR_vis_indices),fitTrackingMatrix(12,ASR_vis_indices),'.','LineStyle','None','MarkerSize',10);
title('ASR Visual \lambda Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,3);
loglog(fitTrackingMatrix(6,ASR_vis_indices),fitTrackingMatrix(13,ASR_vis_indices),'.','LineStyle','None','MarkerSize',10);
title('ASR Visual \omega_0 Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,4);
loglog(fitTrackingMatrix(9,ASR_vis_indices),fitTrackingMatrix(16,ASR_vis_indices),'.','LineStyle','None','MarkerSize',10);
title('ASR Visual f Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

% ASR Auditory
figure;
subplot(2,2,1);
loglog(fitTrackingMatrix(4,ASR_aud_indices),fitTrackingMatrix(11,ASR_aud_indices),'.','LineStyle','None','MarkerSize',10);
title('ASR Auditory A_0 Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,2);
loglog(fitTrackingMatrix(5,ASR_aud_indices),fitTrackingMatrix(12,ASR_aud_indices),'.','LineStyle','None','MarkerSize',10);
title('ASR Auditory \lambda Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,3);
loglog(fitTrackingMatrix(6,ASR_aud_indices),fitTrackingMatrix(13,ASR_aud_indices),'.','LineStyle','None','MarkerSize',10);
title('ASR Auditory \omega_0 Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,4);
loglog(fitTrackingMatrix(9,ASR_aud_indices),fitTrackingMatrix(16,ASR_aud_indices),'.','LineStyle','None','MarkerSize',10);
title('ASR Auditory f Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

% VSR visual
figure;
subplot(2,2,1);
loglog(fitTrackingMatrix(4,VSR_vis_indices),fitTrackingMatrix(11,VSR_vis_indices),'.','LineStyle','None','MarkerSize',10);
title('VSR Visual A_0 Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,2);
loglog(fitTrackingMatrix(5,VSR_vis_indices),fitTrackingMatrix(12,VSR_vis_indices),'.','LineStyle','None','MarkerSize',10);
title('VSR Visual \lambda Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,3);
loglog(fitTrackingMatrix(6,VSR_vis_indices),fitTrackingMatrix(13,VSR_vis_indices),'.','LineStyle','None','MarkerSize',10);
title('VSR Visual \omega_0 Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,4);
loglog(fitTrackingMatrix(9,VSR_vis_indices),fitTrackingMatrix(16,VSR_vis_indices),'.','LineStyle','None','MarkerSize',10);
title('VSR Visual f Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

% VSR Auditory
figure;
subplot(2,2,1);
loglog(fitTrackingMatrix(4,VSR_aud_indices),fitTrackingMatrix(11,VSR_aud_indices),'.','LineStyle','None','MarkerSize',10);
title('VSR Auditory A_0 Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,2);
loglog(fitTrackingMatrix(5,VSR_aud_indices),fitTrackingMatrix(12,VSR_aud_indices),'.','LineStyle','None','MarkerSize',10);
title('VSR Auditory \lambda Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,3);
loglog(fitTrackingMatrix(6,VSR_aud_indices),fitTrackingMatrix(13,VSR_aud_indices),'.','LineStyle','None','MarkerSize',10);
title('VSR Auditory \omega_0 Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,4);
loglog(fitTrackingMatrix(9,VSR_aud_indices),fitTrackingMatrix(16,VSR_aud_indices),'.','LineStyle','None','MarkerSize',10);
title('VSR Auditory f Fitted Values vs Guesses','FontSize',20);
xlabel('Guess','FontSize',16);
ylabel('Fitted Value','FontSize',16);
set(gca,'FontSize',14);

%% ASR Visual ----------------------------------------
fprintf('\n\nASR Visual: \n\n');
disp(ASR_vis);
fprintf('\n\nBest guesses: \n\n');
disp(ASR_vis_guesses);
ASR_vis_A0 = ASR_vis(1,:);
ASR_vis_lambda = ASR_vis(2,:);
ASR_vis_w0 = ASR_vis(3,:);
ASR_vis_f = ASR_vis(6,:);
ASR_vis_B = ASR_vis(7,:);
ASR_vis_A0_w0 = ASR_vis_A0./ASR_vis_w0;

ASR_vis_guess_A0 = ASR_vis_guesses(1,:);
ASR_vis_guess_lambda = ASR_vis_guesses(2,:);
ASR_vis_guess_w0 = ASR_vis_guesses(3,:);
ASR_vis_guess_f = ASR_vis_guesses(6,:);
ASR_vis_guess_B = ASR_vis_guesses(7,:);

ASR_vis_closeParams_A0 = ASR_vis_closeParams(1,:);
ASR_vis_closeParams_lambda = ASR_vis_closeParams(2,:);
ASR_vis_closeParams_w0 = ASR_vis_closeParams(3,:);
ASR_vis_closeParams_f = ASR_vis_closeParams(6,:);
ASR_vis_closeParams_B = ASR_vis_closeParams(7,:);
ASR_vis_closeParams_A0_w0 = ASR_vis_closeParams_A0./ASR_vis_closeParams_w0;

figure;
subplot(2,3,1);
histogram(ASR_vis_guess_A0,5);
title('ASR Visual A_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,2);
histogram(ASR_vis_guess_lambda,5);
title('ASR Visual A_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,3);
title('ASR Visual \omega_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,4);
histogram(ASR_vis_guess_f,5);
title('ASR Visual f guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,5);
histogram(ASR_vis_guess_B,5);
title('ASR Visual B guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

figure;
subplot(2,2,1);
hold on;
histogram(ASR_vis_closeParams_A0_w0,10);
histogram(ASR_vis_A0_w0,5);
title('ASR vis A_0/\omega_0','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

subplot(2,2,2);
hold on;
histogram(ASR_vis_closeParams_lambda,10);
histogram(ASR_vis_lambda,5);
title('ASR vis \lambda','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

subplot(2,2,3);
hold on;
histogram(ASR_vis_closeParams_f,10);
histogram(ASR_vis_f,5);
title('ASR vis f','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

subplot(2,2,4);
hold on;
histogram(ASR_vis_closeParams_B,10);
histogram(ASR_vis_B,5);
title('ASR vis B','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

%% ASR auditory, free S --------------------------------
fprintf('\n\nASR Visual: \n\n');
disp(ASR_aud);
fprintf('\n\nBest guesses: \n\n');
disp(ASR_aud_guesses);
ASR_aud_A0 = ASR_aud(1,:);
ASR_aud_lambda = ASR_aud(2,:);
ASR_aud_w0 = ASR_aud(3,:);
ASR_aud_m = ASR_aud(4,:);
ASR_aud_s = ASR_aud(5,:);
ASR_aud_f = ASR_aud(6,:);
ASR_aud_B = ASR_aud(7,:);
ASR_aud_A0_w0 = ASR_aud_A0./ASR_aud_w0;

ASR_aud_guess_A0 = ASR_aud_guesses(1,:);
ASR_aud_guess_lambda = ASR_aud_guesses(2,:);
ASR_aud_guess_w0 = ASR_aud_guesses(3,:);
ASR_aud_guess_m = ASR_aud_guesses(4,:);
ASR_aud_guess_s = ASR_aud_guesses(5,:)
ASR_aud_guess_f = ASR_aud_guesses(6,:);
ASR_aud_guess_B = ASR_aud_guesses(7,:);

ASR_aud_closeParams_A0 = ASR_aud_closeParams(1,:);
ASR_aud_closeParams_lambda = ASR_aud_closeParams(2,:);
ASR_aud_closeParams_w0 = ASR_aud_closeParams(3,:);
ASR_aud_closeParams_m = ASR_aud_closeParams(4,:);
ASR_aud_closeParams_s = ASR_aud_closeParams(5,:);
ASR_aud_closeParams_f = ASR_aud_closeParams(6,:);
ASR_aud_closeParams_B = ASR_aud_closeParams(7,:);
ASR_aud_closeParams_A0_w0 = ASR_aud_closeParams_A0./ASR_aud_closeParams_w0;

figure;
subplot(2,3,1);
histogram(ASR_aud_guess_A0,5);
title('ASR Visual A_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Guesses within 99% of best','Best Guesses');
set(gca,'FontSize',14);

subplot(2,3,2);
histogram(ASR_aud_guess_lambda,5);
title('ASR Visual A_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,3);
title('ASR Visual \omega_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,4);
histogram(ASR_aud_guess_f,5);
title('ASR Visual f guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,5);
histogram(ASR_aud_guess_B,5);
title('ASR Visual B guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

figure;
subplot(2,3,1);
hold on;
histogram(ASR_aud_closeParams_A0_w0,10);
histogram(ASR_aud_A0_w0,5);
title('ASR aud A_0/\omega_0','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

subplot(2,3,2);
hold on;
histogram(ASR_aud_closeParams_lambda,10);
histogram(ASR_aud_lambda,5);
title('ASR aud \lambda','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

subplot(2,3,3);
hold on;
histogram(ASR_aud_closeParams_m,10);
histogram(ASR_aud_m,5);
title('ASR aud m','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

subplot(2,3,4);
hold on;
histogram(ASR_aud_closeParams_s,10);
histogram(ASR_aud_s,5);
title('ASR aud s','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

subplot(2,3,5);
hold on;
histogram(ASR_aud_closeParams_f,10);
histogram(ASR_aud_f,5);
title('ASR aud f','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);

subplot(2,3,6);
hold on;
histogram(ASR_aud_closeParams_B,10);
histogram(ASR_aud_B,5);
title('ASR aud B','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Within 99% of best','Best');
set(gca,'FontSize',14);


%% ASR tactile ------------------------------------------
% fprintf('\n\nASR Visual: \n\n');
% disp(ASR_tac);
% fprintf('\n\nBest guesses: \n\n');
% disp(ASR_tac_guesses);
% ASR_tac_A0 = ASR_tac(1,:);
% ASR_tac_lambda = ASR_tac(2,:);
% ASR_tac_w0 = ASR_tac(3,:);
% ASR_tac_f = ASR_tac(6,:);
% ASR_tac_B = ASR_tac(7,:);
% ASR_tac_A0_w0 = ASR_tac_A0./ASR_tac_w0;
%
% ASR_tac_guess_A0 = ASR_tac_guesses(1,:);
% ASR_tac_guess_lambda = ASR_tac_guesses(2,:);
% ASR_tac_guess_w0 = ASR_tac_guesses(3,:);
% ASR_tac_guess_f = ASR_tac_guesses(6,:);
% ASR_tac_guess_B = ASR_tac_guesses(7,:);
%
% ASR_tac_closeParams_A0 = ASR_tac_closeParams(1,:);
% ASR_tac_closeParams_lambda = ASR_tac_closeParams(2,:);
% ASR_tac_closeParams_w0 = ASR_tac_closeParams(3,:);
% ASR_tac_closeParams_f = ASR_tac_closeParams(6,:);
% ASR_tac_closeParams_B = ASR_tac_closeParams(7,:);
% ASR_tac_closeParams_A0_w0 = ASR_tac_closeParams_A0./ASR_tac_closeParams_w0;
%
% figure;
% subplot(2,3,1);
% histogram(ASR_tac_guess_A0,5);
% title('ASR Visual A_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% legend('Guesses within 99% of best','Best Guesses');
% set(gca,'FontSize',14);
%
% subplot(2,3,2);
% histogram(ASR_tac_guess_lambda,5);
% title('ASR Visual A_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,3);
% title('ASR Visual \omega_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,4);
% histogram(ASR_tac_guess_f,5);
% title('ASR Visual f guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,5);
% histogram(ASR_tac_guess_B,5);
% title('ASR Visual B guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% figure;
% subplot(2,2,1);
% hold on;
% histogram(ASR_tac_closeParams_A0_w0,10);
% histogram(ASR_tac_A0_w0,5);
% title('ASR tac A_0/\omega_0','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,2);
% hold on;
% histogram(ASR_tac_closeParams_lambda,10);
% histogram(ASR_tac_lambda,5);
% title('ASR tac \lambda','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,3);
% hold on;
% histogram(ASR_tac_closeParams_f,10);
% histogram(ASR_tac_f,5);
% title('ASR tac f','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,4);
% hold on;
% histogram(ASR_tac_closeParams_B,10);
% histogram(ASR_tac_B,5);
% title('ASR tac B','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);

%% ASR vestibular -------------------------------------
% fprintf('\n\nASR Visual: \n\n');
% disp(ASR_vst);
% fprintf('\n\nBest guesses: \n\n');
% disp(ASR_vst_guesses);
% ASR_vst_A0 = ASR_vst(1,:);
% ASR_vst_lambda = ASR_vst(2,:);
% ASR_vst_w0 = ASR_vst(3,:);
% ASR_vst_f = ASR_vst(6,:);
% ASR_vst_B = ASR_vst(7,:);
% ASR_vst_A0_w0 = ASR_vst_A0./ASR_vst_w0;
%
% ASR_vst_guess_A0 = ASR_vst_guesses(1,:);
% ASR_vst_guess_lambda = ASR_vst_guesses(2,:);
% ASR_vst_guess_w0 = ASR_vst_guesses(3,:);
% ASR_vst_guess_f = ASR_vst_guesses(6,:);
% ASR_vst_guess_B = ASR_vst_guesses(7,:);
%
% ASR_vst_closeParams_A0 = ASR_vst_closeParams(1,:);
% ASR_vst_closeParams_lambda = ASR_vst_closeParams(2,:);
% ASR_vst_closeParams_w0 = ASR_vst_closeParams(3,:);
% ASR_vst_closeParams_f = ASR_vst_closeParams(6,:);
% ASR_vst_closeParams_B = ASR_vst_closeParams(7,:);
% ASR_vst_closeParams_A0_w0 = ASR_vst_closeParams_A0./ASR_vst_closeParams_w0;
%
% figure;
% subplot(2,3,1);
% histogram(ASR_vst_guess_A0,5);
% title('ASR Visual A_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% legend('Guesses within 99% of best','Best Guesses');
% set(gca,'FontSize',14);
%
% subplot(2,3,2);
% histogram(ASR_vst_guess_lambda,5);
% title('ASR Visual A_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,3);
% title('ASR Visual \omega_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,4);
% histogram(ASR_vst_guess_f,5);
% title('ASR Visual f guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,5);
% histogram(ASR_vst_guess_B,5);
% title('ASR Visual B guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% figure;
% subplot(2,2,1);
% hold on;
% histogram(ASR_vst_closeParams_A0_w0,10);
% histogram(ASR_vst_A0_w0,5);
% title('ASR vst A_0/\omega_0','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,2);
% hold on;
% histogram(ASR_vst_closeParams_lambda,10);
% histogram(ASR_vst_lambda,5);
% title('ASR vst \lambda','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,3);
% hold on;
% histogram(ASR_vst_closeParams_f,10);
% histogram(ASR_vst_f,5);
% title('ASR vst f','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,4);
% hold on;
% histogram(ASR_vst_closeParams_B,10);
% histogram(ASR_vst_B,5);
% title('ASR vst B','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);

%% VSR visual ----------------------------------------
fprintf('\n\nASR Visual: \n\n');
disp(VSR_vis);
fprintf('\n\nBest guesses: \n\n');
disp(VSR_vis_guesses);
VSR_vis_A0 = VSR_vis(1,:);
VSR_vis_lambda = VSR_vis(2,:);
VSR_vis_w0 = VSR_vis(3,:);
VSR_vis_f = VSR_vis(6,:);
VSR_vis_B = VSR_vis(7,:);
VSR_vis_A0_w0 = VSR_vis_A0./VSR_vis_w0;

VSR_vis_guess_A0 = VSR_vis_guesses(1,:);
VSR_vis_guess_lambda = VSR_vis_guesses(2,:);
VSR_vis_guess_w0 = VSR_vis_guesses(3,:);
VSR_vis_guess_f = VSR_vis_guesses(6,:);
VSR_vis_guess_B = VSR_vis_guesses(7,:);

VSR_vis_closeParams_A0 = VSR_vis_closeParams(1,:);
VSR_vis_closeParams_lambda = VSR_vis_closeParams(2,:);
VSR_vis_closeParams_w0 = VSR_vis_closeParams(3,:);
VSR_vis_closeParams_f = VSR_vis_closeParams(6,:);
VSR_vis_closeParams_B = VSR_vis_closeParams(7,:);
VSR_vis_closeParams_A0_w0 = VSR_vis_closeParams_A0./VSR_vis_closeParams_w0;

figure;
subplot(2,3,1);
histogram(VSR_vis_guess_A0,5);
title('ASR Visual A_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Guesses within 99% of best','Best Guesses');
set(gca,'FontSize',14);

subplot(2,3,2);
histogram(VSR_vis_guess_lambda,5);
title('ASR Visual A_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,3);
title('ASR Visual \omega_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,4);
histogram(VSR_vis_guess_f,5);
title('ASR Visual f guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,5);
histogram(VSR_vis_guess_B,5);
title('ASR Visual B guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

figure;
subplot(2,2,1);
hold on;
histogram(VSR_vis_closeParams_A0_w0,10);
histogram(VSR_vis_A0_w0,5);
title('VSR vis A_0/\omega_0','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,2);
hold on;
histogram(VSR_vis_closeParams_lambda,10);
histogram(VSR_vis_lambda,5);
title('VSR vis \lambda','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,3);
hold on;
histogram(VSR_vis_closeParams_f,10);
histogram(VSR_vis_f,5);
title('VSR vis f','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,4);
hold on;
histogram(VSR_vis_closeParams_B,10);
histogram(VSR_vis_B,5);
title('VSR vis B','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

%% VSR Auditory ---------------------------------------
fprintf('\n\nASR Visual: \n\n');
disp(VSR_aud);
fprintf('\n\nBest guesses: \n\n');
disp(VSR_aud_guesses);
VSR_aud_A0 = VSR_aud(1,:);
VSR_aud_lambda = VSR_aud(2,:);
VSR_aud_w0 = VSR_aud(3,:);
VSR_aud_f = VSR_aud(6,:);
VSR_aud_B = VSR_aud(7,:);
VSR_aud_A0_w0 = VSR_aud_A0./VSR_aud_w0;

VSR_aud_guess_A0 = VSR_aud_guesses(1,:);
VSR_aud_guess_lambda = VSR_aud_guesses(2,:);
VSR_aud_guess_w0 = VSR_aud_guesses(3,:);
VSR_aud_guess_f = VSR_aud_guesses(6,:);
VSR_aud_guess_B = VSR_aud_guesses(7,:);

VSR_aud_closeParams_A0 = VSR_aud_closeParams(1,:);
VSR_aud_closeParams_lambda = VSR_aud_closeParams(2,:);
VSR_aud_closeParams_w0 = VSR_aud_closeParams(3,:);
VSR_aud_closeParams_f = VSR_aud_closeParams(6,:);
VSR_aud_closeParams_B = VSR_aud_closeParams(7,:);
VSR_aud_closeParams_A0_w0 = VSR_aud_closeParams_A0./VSR_aud_closeParams_w0;

figure;
subplot(2,3,1);
histogram(VSR_aud_guess_A0,5);
title('ASR Visual A_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
legend('Guesses within 99% of best','Best Guesses');
set(gca,'FontSize',14);

subplot(2,3,2);
histogram(VSR_aud_guess_lambda,5);
title('ASR Visual A_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,3);
title('ASR Visual \omega_0 guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,4);
histogram(VSR_aud_guess_f,5);
title('ASR Visual f guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,3,5);
histogram(VSR_aud_guess_B,5);
title('ASR Visual B guesses','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

figure;
subplot(2,2,1);
hold on;
histogram(VSR_aud_closeParams_A0_w0,10);
histogram(VSR_aud_A0_w0,5);
title('VSR aud A_0/\omega_0','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,2);
hold on;
histogram(VSR_aud_closeParams_lambda,10);
histogram(VSR_aud_lambda,5);
title('VSR aud \lambda','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,3);
hold on;
histogram(VSR_aud_closeParams_f,10);
histogram(VSR_aud_f,5);
title('VSR aud f','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

subplot(2,2,4);
hold on;
histogram(VSR_aud_closeParams_B,10);
histogram(VSR_aud_B,5);
title('VSR aud B','FontSize',20);
xlabel('Value','FontSize',16);
ylabel('Number of Occurrences','FontSize',16);
set(gca,'FontSize',14);

%% VSR tactile ------------------------------------
% fprintf('\n\nASR Visual: \n\n');
% disp(VSR_tac);
% fprintf('\n\nBest guesses: \n\n');
% disp(VSR_tac_guesses);
% VSR_tac_A0 = VSR_tac(1,:);
% VSR_tac_lambda = VSR_tac(2,:);
% VSR_tac_w0 = VSR_tac(3,:);
% VSR_tac_f = VSR_tac(6,:);
% VSR_tac_B = VSR_tac(7,:);
% VSR_tac_A0_w0 = VSR_tac_A0./VSR_tac_w0;
%
% VSR_tac_guess_A0 = VSR_tac_guesses(1,:);
% VSR_tac_guess_lambda = VSR_tac_guesses(2,:);
% VSR_tac_guess_w0 = VSR_tac_guesses(3,:);
% VSR_tac_guess_f = VSR_tac_guesses(6,:);
% VSR_tac_guess_B = VSR_tac_guesses(7,:);
%
% VSR_tac_closeParams_A0 = VSR_tac_closeParams(1,:);
% VSR_tac_closeParams_lambda = VSR_tac_closeParams(2,:);
% VSR_tac_closeParams_w0 = VSR_tac_closeParams(3,:);
% VSR_tac_closeParams_f = VSR_tac_closeParams(6,:);
% VSR_tac_closeParams_B = VSR_tac_closeParams(7,:);
% VSR_tac_closeParams_A0_w0 = VSR_tac_closeParams_A0./VSR_tac_closeParams_w0;
%
% figure;
% subplot(2,3,1);
% histogram(VSR_tac_guess_A0,5);
% title('ASR Visual A_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% legend('Guesses within 99% of best','Best Guesses');
% set(gca,'FontSize',14);
%
% subplot(2,3,2);
% histogram(VSR_tac_guess_lambda,5);
% title('ASR Visual A_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,3);
% title('ASR Visual \omega_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,4);
% histogram(VSR_tac_guess_f,5);
% title('ASR Visual f guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,5);
% histogram(VSR_tac_guess_B,5);
% title('ASR Visual B guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% figure;
% subplot(2,2,1);
% hold on;
% histogram(VSR_tac_closeParams_A0_w0,10);
% histogram(VSR_tac_A0_w0,5);
% title('VSR tac A_0/\omega_0','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,2);
% hold on;
% histogram(VSR_tac_closeParams_lambda,10);
% histogram(VSR_tac_lambda,5);
% title('VSR tac \lambda','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,3);
% hold on;
% histogram(VSR_tac_closeParams_f,10);
% histogram(VSR_tac_f,5);
% title('VSR tac f','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,4);
% hold on;
% histogram(VSR_tac_closeParams_B,10);
% histogram(VSR_tac_B,5);
% title('VSR tac B','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);

%% VSR vestibular ----------------------------------------
% fprintf('\n\nASR Visual: \n\n');
% disp(VSR_vst);
% fprintf('\n\nBest guesses: \n\n');
% disp(VSR_vst_guesses);
% VSR_vst_A0 = VSR_vst(1,:);
% VSR_vst_lambda = VSR_vst(2,:);
% VSR_vst_w0 = VSR_vst(3,:);
% VSR_vst_f = VSR_vst(6,:);
% VSR_vst_B = VSR_vst(7,:);
% VSR_vst_A0_w0 = VSR_vst_A0./VSR_vst_w0;
%
% VSR_vst_guess_A0 = VSR_vst_guesses(1,:);
% VSR_vst_guess_lambda = VSR_vst_guesses(2,:);
% VSR_vst_guess_w0 = VSR_vst_guesses(3,:);
% VSR_vst_guess_f = VSR_vst_guesses(6,:);
% VSR_vst_guess_B = VSR_vst_guesses(7,:);
%
% VSR_vst_closeParams_A0 = VSR_vst_closeParams(1,:);
% VSR_vst_closeParams_lambda = VSR_vst_closeParams(2,:);
% VSR_vst_closeParams_w0 = VSR_vst_closeParams(3,:);
% VSR_vst_closeParams_f = VSR_vst_closeParams(6,:);
% VSR_vst_closeParams_B = VSR_vst_closeParams(7,:);
% VSR_vst_closeParams_A0_w0 = VSR_vst_closeParams_A0./VSR_vst_closeParams_w0;
%
% figure;
% subplot(2,3,1);
% histogram(VSR_vst_guess_A0,5);
% title('ASR Visual A_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% legend('Guesses within 99% of best','Best Guesses');
% set(gca,'FontSize',14);
%
% subplot(2,3,2);
% histogram(VSR_vst_guess_lambda,5);
% title('ASR Visual A_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,3);
% title('ASR Visual \omega_0 guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,4);
% histogram(VSR_vst_guess_f,5);
% title('ASR Visual f guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,3,5);
% histogram(VSR_vst_guess_B,5);
% title('ASR Visual B guesses','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% figure;
% subplot(2,2,1);
% hold on;
% histogram(VSR_vst_closeParams_A0_w0,10);
% histogram(VSR_vst_A0_w0,5);
% title('VSR vst A_0/\omega_0','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,2);
% hold on;
% histogram(VSR_vst_closeParams_lambda,10);
% histogram(VSR_vst_lambda,5);
% title('VSR vst \lambda','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,3);
% hold on;
% histogram(VSR_vst_closeParams_f,10);
% histogram(VSR_vst_f,5);
% title('VSR vst f','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
%
% subplot(2,2,4);
% hold on;
% histogram(VSR_vst_closeParams_B,10);
% histogram(VSR_vst_B,5);
% title('VSR vst B','FontSize',20);
% xlabel('Value','FontSize',16);
% ylabel('Number of Occurrences','FontSize',16);
% set(gca,'FontSize',14);
% ASR Aud S fits ----------------------------------------

figure;

subplot(2,1,1);
histogram(ASR_aud_s,5);
title('ASR Auditory Start of Masking (s)','FontSize',20);
axis([10 30 0 5]);
set(gca,'FontSize',14);

subplot(2,1,2);
histogram(ASR_aud_s - ASR_aud_B,5);
title('Difference Between Start of Masking and Baseline');

runtime = toc;
fprintf('Program Run Time: %4.3f minutes\n',runtime/60);