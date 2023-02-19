function [allUnderlyingThresholds,allMeasuredThresholds,A0Guess,lambdaGuess,w0Guess,mGuess,sGuess,fGuess,BGuess,B_sd,Ymin,Ymax,underlyingCurveParameters ,UpperBounds,LowerBounds,options,...
    guess6,BTL,MST] = curveFitSetup(SRType,TestType,maskingExpected,levels,thisSubjectMeasuredThresholds)


%% Globals
global levels thisSubjectMeasuredThresholds TestType SRType method

%% Task Setup

% Best Threshold Location
BTL = levels(1+min(find(thisSubjectMeasuredThresholds(2:end) == min(thisSubjectMeasuredThresholds(2:end)))))+0.00001;

% Measured Sham Threshold
MST = thisSubjectMeasuredThresholds(1);

% Matrix of all underlying parameters
%allUnderlyingParameters = zeros(7,numberOfSimulations);

% Matrix of all fitted parameters
%allFittedParameters = zeros(7,numberOfSimulations);

% Matrix of underlying thresholds
allUnderlyingThresholds = [];

% Matrix of fitted thresholds
allMeasuredThresholds = [];

% The feature names cell array
featureNames = {};

% The underlying parameters for subjects who exhibit SR. May be changed for
% different task types.
switch SRType
    
    case 1
        
        switch TestType
            
            case 1
                
                % ASR Visual
                A0Guess = [0.342;0;1;10;1];
                lambdaGuess = [5;0;10;1;5];
                w0Guess = [5;0;10;1;10];
                mGuess = [0;0;0;0;0];
                sGuess = [0;0;0;0;0];
                fGuess = [BTL;BTL;BTL;BTL;BTL];
                f_vec = (5.01:1:15.01);
                BGuess = [MST;MST;MST;MST;MST];
                
                % The expected standard deviation of baseline thresholds
                B_sd = 0.05;
                
                % Curve plot axes limits
                Ymin = 0;
                Ymax = 1;
                
                % Curve Fit Guesses
                guess6 = [0;0;0;0;0;0;0];
                
            case 2
                
      % ASR Auditory
                A0Guess = [244.2;0;10;100;100];
                lambdaGuess = [3;0;10;5;50];
                w0Guess = [50;0;10;1;1];
                m = [0.6;0.6;0.6;0.6;0.6];
                SGuess = [(MST+15);MST+15;MST+15;MST+15;MST+15];               
                fGuess = [(MST-5);MST-5;MST-5;MST-5;MST-5];
                f_vec = (5.01:1:15.01);
                BGuess = [MST;MST;MST;MST;MST];
                
                % Expected standard deviation of baseline thresholds
                B_sd = 1;
                
                % Curve plot axes
                Ymin = 0;
                Ymax = 20;
                
                % Curve Fit Guesses
                guess6 = [0;0;0;0.6;MST+15;MST-5;MST];
                
            case 3
                
                % ASR Tactile
                A0Guess = [0.615;0;1;10;1];
                lambdaGuess = [5;0;10;1;5];
                w0Guess = [5;0;10;1;10];
                mGuess = [0;0;0;0;0];
                sGuess = [0;0;0;0;0];              
                fGuess = [BTL;BTL;BTL;BTL;BTL];
                f_vec = (5.01:1:15.01);
                BGuess = [MST;MST;MST;MST;MST];
                
                % Expected standard deviation of baseline thresholds
                B_sd = 0.1;
                
                % Curve plot axes
                Ymin = 0;
                Ymax = 1;
                
                % Curve Fit Guesses
                guess6 = [0;0;0;0;0;0;0];
                
        end
        
    case 2 % VSR
        
        switch TestType
            
            case 1
                
                % VSR Visual
                A0Guess = [0.342;0;0.2;0.1;0.5];
                lambdaGuess = [0.1;0;0.2;0.3;0.05];
                w0Guess = [0.1;0;0.2;0.3;0.05];
                mGuess = [0;0;0;0;0];
                sGuess = [0;0;0;0;0];
                fGuess = [BTL;BTL;BTL;BTL;BTL];
                f_vec = [0.201:0.1:0.601];
                BGuess = [MST;MST;MST;MST;MST];
                
                % Expected standard deviation of baseline thresholds
                B_sd = 0.05;
                
                % Curve plot axes
                Ymin = 0;
                Ymax = 0.3;
                
                % Curve Fit Guesses
                guess6 = [0;0;0;0;0;0;0];
                
            case 2
                
                % VSR Auditory
                A0Guess = [15.158;0;10;5;25];
                lambdaGuess = [0.1;0;0.2;0.3;0.05];
                w0Guess = [0.1;0;0.2;0.3;0.05];
                mGuess = [0;0;0;0;0];
                sGuess = [0;0;0;0;0];
                fGuess = [BTL;BTL;BTL;BTL;BTL];
                f_vec = [0.201:0.1:0.601];
                BGuess = [MST;MST;MST;MST;MST];
                
                % Expected standard deviation of baseline thresholds
                B_sd = 1.5;
                
                % Curve plot axes
                Ymin = -5;
                Ymax = 20;
                
                % Curve Fit Guesses
                guess6 = [0;0;0;0;0;0;0];
                
                
            case 3
                
                % VSR Tactile
                A0Guess = [0.615;0;10;5;25];
                lambdaGuess = [0.1;0;0.2;0.3;0.05];
                w0Guess = [0.1;0;0.2;0.3;0.05];
                mGuess = [0;0;0;0;0];
                sGuess = [0;0;0;0;0];
                fGuess = [BTL;BTL;BTL;BTL;BTL];
                f_vec = [0.201:0.1:0.601];
                BGuess = [MST;MST;MST;MST;MST];
                
                % Expected standard deviation of baseline thresholds
                B_sd = 0.1;
                
                % Curve plot axes
                Ymin = 0;
                Ymax = 1;
                
                % Curve Fit Guesses
                guess6 = [0;0;0;0;0;0;0];
        end
end
% Generate the underlying curve parameters
underlyingCurveParameters = [A0Guess;
    lambdaGuess;
    w0Guess;
    mGuess;
    sGuess;
    fGuess;
    BGuess];

%% Guesses for curve fits
if maskingExpected % Case with masking
    
    % Use all the parameters. Must be modified slightly so initial guesses
    % don't cause infinite values in first evaluation
    curveFitGuesses = underlyingCurveParameters;
    
    
else % No masking
    
    % Don't use parameters that define masking. Must be modified slightly
    % so initial guesses don't cause infinite values in first evaluation.
    curveFitGuesses = [underlyingCurveParameters(1:3);underlyingCurveParameters(6:7)];
    
end

%% Upper and lower bounds on curve fit parameters
if maskingExpected % Case with masking
    
    % Upper bounds on all 7 parameters
    UpperBounds = [1e10;
        1e10;
        1e10;
        1e10;
        1e10;
        1e10;
        1e10];
    
    % Lower bounds on all 7 parameters
    LowerBounds = [1e-6;
        1e-6;
        1e-6;
        1e-6;
        1e-6;
        1e-6;
        1e-6];
    
else % No masking
    
    % Upper bounds on 5 relevant parameters
    UpperBounds = [1e10;
        1e10;
        1e10;
        1e10;
        1e10];
    
    % Lower bounds on 5 relevant parameters
    LowerBounds = [1e-6;
        1e-6;
        1e-6;
        1e-6;
        1e-6];
end

%% Curve fit options

% [USER EDIT]: Change the display options to show whether curve fits
% converge
options = optimoptions('fmincon','Display','off');
% [USER EDIT]: Change warning display
warning('off');


% ===== END SETUP | BEGIN SIMULATION ====================================================================================================================
