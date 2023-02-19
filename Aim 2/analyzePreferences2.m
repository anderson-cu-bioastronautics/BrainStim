function [mdlD,NoisePreferN_mat,testDataN_mat] = analyzePreferences2(NoisePrefer_mat,testData_mat);

%% Regress the data and define the new data vectors

% Feedback data
for kay = 1:3
    NoisePrefer = NoisePrefer_mat{kay}; testData{kay} = testData_mat{kay};
    mdl = fitlm(NoisePrefer, testData{kay}, 'linear');
    t_cookd = 3*mean(mdl.Diagnostics.CooksDistance,'omitnan');
    cookBad = find(mdl.Diagnostics.CooksDistance > t_cookd);
    testBad{1,kay} = cookBad;
    mdlD{1,kay} = mdl;
    mdl=[]; cookbad=[];
end

%% Define new vectors based on outlier analysis

for kay=1:3
    % Create new data matrices
    NoisePreferN = NoisePrefer_mat{kay};
    testDataN = testData{kay};
    NoisePreferN(testBad{kay}) = [];
    testDataN(testBad{kay}) = [];
    NoisePreferN_mat{1,kay} = NoisePreferN;
    testDataN_mat{1,kay} = testDataN;
    NoisePreferN = []; testDataN = [];
        
end