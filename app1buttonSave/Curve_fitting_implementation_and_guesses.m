function [curveFitParameters,A0_fit,lambda_fit,w0_fit,m_fit,s_fit,f_fit,B_fit,X,r,SRDip,mask,SRCurveFit]...
    = Curve_fitting_implementation_and_guesses(underlyingCurveParameters,maskingExpected,levels,thisSubjectMeasuredThresholds,LowerBounds,UpperBounds,...
    options,underlyingShamThreshold,A0Guess,lambdaGuess,w0Guess,mGuess,sGuess,fGuess,BGuess,B_sd,guess6,BTL,MST);

%% Curve Fitting
    
global TestType SRType levels thisSubjectMeasuredThresholds method widthThreshold1 widthThreshold2

    % Re-establish the cost function to be really high
    Jmin = 1e10;
    
    thisSubjectSE = ones(size(levels));
    
% Determine whether masking is expected
if maskingExpected
       
    %% Curve fits - With Masking - SR Allowed

%     for guessLoop = 1:length(A0Guess)
      % Guess 1   
        A0 = A0Guess(1);
        lambda = lambdaGuess(1);
        w0 = w0Guess(1);
        m = mGuess(1);
        s = sGuess(1);
        f = fGuess(1);
        B = BGuess(1);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_withmask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
            
    % Guess 2
        A0 = A0Guess(2);
        lambda = lambdaGuess(2);
        w0 = w0Guess(2);
        m = mGuess(2);
        s = sGuess(2);
        f = fGuess(2);
        B = BGuess(2);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_withmask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
       
       % Guess 3
        A0 = A0Guess(3);
        lambda = lambdaGuess(3);
        w0 = w0Guess(3);
        m = mGuess(3);
        s = sGuess(3);
        f = fGuess(3);
        B = BGuess(3);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_withmask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
                 
       % Guess 4
        A0 = A0Guess(4);
        lambda = lambdaGuess(4);
        w0 = w0Guess(4);
        m = mGuess(4);
        s = sGuess(4);
        f = fGuess(4);
        B = BGuess(4);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_withmask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
            
       % Guess 5
        A0 = A0Guess(5);
        lambda = lambdaGuess(5);
        w0 = w0Guess(5);
        m = mGuess(5);
        s = sGuess(5);
        f = fGuess(5);
        B = BGuess(5);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_withmask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
                 
%     end

    %% Curve fits - With Masking - No SR Allowed - Guess 6
       
          [curveFitRawParameters_noSR_allowed,Jnew,flag_new] = fminsearch(@(curveFitRawParameters_noSR_allowed) SR_curve_fit_noSR_withmask(curveFitRawParameters_noSR_allowed,thisSubejctSRNoiseLevels,thisSubjectMeasuredThresholds,thisSubjectSE),guess6);
         
        % Evaluate whether this curve fit is better
        if Jnew < Jmin
            
            % Assign the new cost function value
            Jmin = Jnew;
            
            curveFitRawParameters = curveFitRawParameters_noSR_allowed;
            
            % Extract all the fitted values
            A0_fit = 0;
            lambda_fit = 1;
            w0_fit = 1;
            m_fit = curveFitRawParameters(1);
            s_fit = curveFitRawParameters(2);
            f_fit = 0.000001;
            B_fit = curveFitRawParameters(3);
            
            % Cost function value for printing
            J = Jnew;
            
            % Exit flag value for printing
            flag = flag_new;
           
        end
        
else
        
    %% Do the curve fits - No Masking - SR Allowed

    % Guess 1   
        A0 = A0Guess(1);
        lambda = lambdaGuess(1);
        w0 = w0Guess(1);
        m = mGuess(1);
        s = sGuess(1);
        f = fGuess(1);
        B = BGuess(1);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_nomask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
            
    % Guess 2
        A0 = A0Guess(2);
        lambda = lambdaGuess(2);
        w0 = w0Guess(2);
        m = mGuess(2);
        s = sGuess(2);
        f = fGuess(2);
        B = BGuess(2);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_nomask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
       
       % Guess 3
        A0 = A0Guess(3);
        lambda = lambdaGuess(3);
        w0 = w0Guess(3);
        m = mGuess(3);
        s = sGuess(3);
        f = fGuess(3);
        B = BGuess(3);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_nomask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
                 
       % Guess 4
        A0 = A0Guess(4);
        lambda = lambdaGuess(4);
        w0 = w0Guess(4);
        m = mGuess(4);
        s = sGuess(4);
        f = fGuess(4);
        B = BGuess(4);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_nomask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
            
       % Guess 5
        A0 = A0Guess(5);
        lambda = lambdaGuess(5);
        w0 = w0Guess(5);
        m = mGuess(5);
        s = sGuess(5);
        f = fGuess(5);
        B = BGuess(5);
            
        guesses = [A0;lambda;w0;m;s;f;B];
            
        [curveFitRawParameters_SR,Jnew,flag_new] = fmincon(@(curveFitRawParameters_SR) SR_curve_fit_nomask(curveFitRawParameters_SR,levels,thisSubjectMeasuredThresholds,thisSubjectSE),guesses,[],[],[],[],LowerBounds,UpperBounds,@SR_minWidth_maxDepth_Updated,options);
        % Determine whether this curve fit was better
            if Jnew <= Jmin
                
                % Assign new minimum cost function vlaue
                    Jmin = Jnew;
                % Assign the raw parameters
                    curveFitRawParameters = curveFitRawParameters_SR;
                % Extract all the fitted values
                    A0_fit = curveFitRawParameters(1);
                    lambda_fit = curveFitRawParameters(2);
                    w0_fit = curveFitRawParameters(3);
                    m_fit = 0;
                    s_fit = 0;
                    f_fit = curveFitRawParameters(4);
                    B_fit = curveFitRawParameters(5);
                    
                % Cost function value for printing
                    J = Jnew;
                % Exit flag value for printing
                    flag = flag_new;        
            
            end
         
     %% Do the curve fits - No Masking - No SR allowed
         
        B_fit = sum(thisSubjectMeasuredThresholds./(thisSubjectSE.^2))./length(thisSubjectMeasuredThresholds);
        Jnew = (thisSubjectMeasuredThresholds - B_fit).^2;
        Jnew = sum(Jnew);
        
        % Evaluate whether this curve fit is better
        if Jnew <= Jmin
            
            % Assign new minimum cost function vlaue
            Jmin = Jnew;
            
            % Assign the raw parameters
            
            
            % Extract all the fitted values
            A0_fit = 0;
            lambda_fit = 1;
            w0_fit = 1;
            m_fit = 0;
            s_fit = 0;
            f_fit = 0.00001;
            B_fit = B_fit;
            
            % The cost function value for printing
            J = Jnew;
            
            % Exit flag value for printing
            flag = 1;
            %exitFlags(simulation) = 1;
            
        end
            
end
    
   %% Organize all of the curve fit parameters into a nice, neat vector
    curveFitParameters = [A0_fit;
        lambda_fit;
        w0_fit;
        m_fit;
        s_fit;
        f_fit;
        B_fit];
 
   %% Do the calculations
    X = linspace(min(levels),max(levels),1000);
    
    % Interim calculation of r
    r = (lambda_fit./(sqrt(2)*pi)).*exp(-lambda_fit^2./(2*((X-f_fit).^2)));
    
    % Calculate the SR dip
    SRDip = B_fit - (A0_fit.*lambda_fit./(X - f_fit).^2).*r./sqrt(4.*r.^2 + w0_fit^2).*(X-f_fit).*(X >=f_fit);
    
    % Calculate masking effects
    mask = m_fit*(X-s_fit).*(X >= s_fit);
    
    % Determine the fitted thresholds
    SRCurveFit = SRDip + mask;