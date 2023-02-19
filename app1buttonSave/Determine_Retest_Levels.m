 %% Determine the retest noise levels
    
    % Find the curve fit minimum
    curveFitMinimum = min(X(find(SRCurveFit == min(SRCurveFit))));
  
    if SRType == 1 && TestType == 2
        % Find the rounded
    roundedCFIndex = 1+find(abs(levles(2:end-1) - curveFitMinimum) == min(abs(levles(2:end-1) - curveFitMinimum)));
    
    % Find the closest initial swath level
    roundedCFMinimum = levles(roundedCFIndex);
    
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
        
    elseif roundedCFIndex == length(levles)-1
        
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
    
    retestLevels = levles(retestIndices);
    

    else
    % Find the rounded
    roundedCFIndex = 1+find(abs(levles(2:end) - curveFitMinimum) == min(abs(levles(2:end) - curveFitMinimum)));
    
    % Find the closest initial swath level
    roundedCFMinimum = levles(roundedCFIndex);
    
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
        
    elseif roundedCFIndex == length(levles)
        
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
    retestLevels = levles(retestIndices);
    lengthRetestLevels = length(retestLevels)
    fprintf('The Retest Values are %.2f, %.2f, %.2f, %.2f \n',retestLevels(1),retestLevels(2),retestLevels(3),retestLevels(4))
    
    %% Make sure to randomize the retest levels and say how many trials to do at each!
    
    % Visual/tactile tasks: 100 at sham/best, 50 at neighbors
    % Auditory task: 175 at sham/best, 88 at neighbors
    % ASR Auditory: TBD, so use auditory task for now