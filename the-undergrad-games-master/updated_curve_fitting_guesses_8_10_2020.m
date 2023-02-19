%% Curve Fitting

% Reset the "best" cost function value
Jmin = 1e10;
   
    % Use the minimum measured threshold to determine the curve fit dip
    % location guess
    f_guess = thisSubjectSRNoiseLevels(min(find(thisSubjectMeasuredThresholds(2:end) == min(thisSubjectMeasuredThresholds(2:end)))))+0.001;
   
    s_guess = thisSubjectMeasuredThresholds(1)+15;
   
    B_guess = thisSubjectMeasuredThresholds(1);
   
    % Determine whether masking is expected
    if maskingExpected
       
        curveFitGuesses = [thisSubjectUnderlyingCurveParameters(1:4);s_guess;f_guess;B_guess];
       
        % The fmincon function that returns all the curve fit parameters, for
        % cases where masking is expected
        [curveFitRawParameters_SR,J_SR] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_withmask(curveFitRawParameters_SR,thisSubjectSRNoiseLevels,thisSubjectMeasuredThresholds,thisSubjectSE),curveFitGuesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
       
       
        if J_SR < Jmin
           
            Jmin = J_SR;
           
            curveFitRawParameters = curveFitRawParameters_SR;
           
           
        end

          % Update the guesses to assume no SR
        noSRGuesses = [0;0;0;curveFitGuesses(4:7)];  
       
        % Re-run the curve fit assuming no SR
        [curveFitRawParameters_noSR,J_noSR] = fmincon(@(curveFitRawParameters_noSR) SR_curve_fit_nomask(curveFitRawParameters_noSR,thisSubjectSRNoiseLevels,thisSubjectMeasuredThresholds,thisSubjectSE),noSRGuesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
       
           
        if J_noSR < Jmin
           
            Jmin = J_noSR;
           
            curveFitRawParameters = curveFitRawParameters_noSR;
           
    end

% Extract all the fitted values
        A0_fit = curveFitRawParameters(1);
        lambda_fit = curveFitRawParameters(2);
        w0_fit = curveFitRawParameters(3);
        m_fit = curveFitRawParameters(4);
        s_fit = curveFitRawParameters(5);
        f_fit = curveFitRawParameters(6);
        B_fit = curveFitRawParameters(7);

else % Case with no masking
   
        curveFitGuesses = [thisSubjectUnderlyingCurveParameters(1:3);f_guess;B_guess];
       
        % The fmincon function that returns all the curve fit parameters, for
        % cases where no masking is expected
        [curveFitRawParameters_SR,J_SR] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_nomask(curveFitRawParameters_SR,thisSubjectSRNoiseLevels,thisSubjectMeasuredThresholds,thisSubjectSE),curveFitGuesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
       
       
        if J_SR < Jmin
           
            Jmin = J_SR;
           
            curveFitRawParameters = curveFitRawParameters_SR;
           
        end
        % Update initial guesses to assume no SR
        noSRGuesses = [0;0;0;curveFitGuesses(4:5)];
       
        % Re-run the curve fit assuming no SR
        [curveFitRawParameters_noSR,J_noSR,flag_noSR] = fmincon(@(curveFitRawParameters_noSR) SR_curve_fit_nomask(curveFitRawParameters_noSR,thisSubjectSRNoiseLevels,thisSubjectMeasuredThresholds,thisSubjectSE),noSRGuesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);

        if J_noSR < Jmin
           
            Jmin = J_noSR;
           
            curveFitRawParameters = curveFitRawParameters_noSR;
           

           
        end

% Extract all the fitted values
        A0_fit = curveFitRawParameters(1);
        lambda_fit = curveFitRawParameters(2);
        w0_fit = curveFitRawParameters(3);
        m_fit = 0;
        s_fit = 0;
        f_fit = curveFitRawParameters(4);
        B_fit = curveFitRawParameters(5);

end

% Organize all of the curve fit parameters into a nice, neat vector
    curveFitParameters = [A0_fit;
        lambda_fit;
        w0_fit;
        m_fit;
        s_fit;
        f_fit;
        B_fit];
   
    % Clear the raw parameters for the next iteration
    clear curveFitRawParameters;