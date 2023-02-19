clear all;
close all;

TestType = 'Visual';
SRType = 'ASR';

%get files
searchterm0 = ['../SubjectData/*/*' SRType '*' TestType '*0mA_0dB_2.xls'];
searchtermB = ['../SubjectData/*/*' SRType '*' TestType '*_2.xls'];
listing0 = dir(searchterm0);
listingB = dir(searchtermB);

for i =1:length(listingB)
    filenames{i} = [listingB(i).folder '\' listingB(i).name];
end

for i =1:length(listing0)
    filenamesham{i} = [listing0(i).folder '\' listing0(i).name];
end

filenamesretest = setdiff(filenames, filenamesham);

sham_sigmas = zeros(1, length(filenamesham));
sham_thresholds = zeros(1, length(filenamesham));
retest_sigmas = zeros(1, length(filenamesretest));
retest_thresholds =zeros(1, length(filenamesretest));

if strcmp(TestType, 'Visual')
    %visual
    mu_guess = 0.2; sigma_guess = 0.05;
elseif strcmp(TestType, 'Auditory')
    %auditory
    mu_guess = 8; sigma_guess = 3;
elseif strcmp(TestType, 'Tactile')
    %Tactile
    mu_guess = 0.3; sigma_guess = 0.1;
else
    mu_guess = 'something is wrong'
    sigma_guess = 'something is wrong'
end

for i = 1:length(filenamesham)
    filename = filenamesham{i};
    master=importdata(filename);
    if strcmp(TestType, 'Auditory') 
        X = master.data.Sheet1(:,2);
        cor = master.data.Sheet1(:,4);
    else
        X = master.data.Sheet1(:,2);
        X = abs(X);
        cor = master.data.Sheet1(:,4);
    end
    
    x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mu_guess, sigma_guess]);
    
    sham_thresholds(i) = x(1);
    sham_sigmas(i) = x(2);
end

for i = 1:length(filenamesretest)
    filename = filenamesretest{i};
    master=importdata(filename);
    if strcmp(TestType, 'Auditory')
        X = master.data.Sheet1(:,2);
        cor = master.data.Sheet1(:,4);
    else
        X = master.data.Sheet1(:,2);
        X = abs(X);
        cor = master.data.Sheet1(:,4);
    end
    
    x = fminsearch(@(x) two_int_fit_simp(x, X, cor), [mu_guess, sigma_guess]);
    retest_thresholds(i) = x(1);
    retest_sigmas(i) = x(2);
end



if length(filenamesretest) == length(filenamesham)
    [h,p,ci,stats] = ttest(sham_thresholds, retest_thresholds)
else
    [h,p,ci,stats] = ttest2(sham_thresholds,retest_thresholds)
end

[h,p,ci,stats]=ttest(sham_thresholds, retest_thresholds,'Tail','right','Alpha',0.05)

 
x = [0 1];
y = [mean(sham_thresholds) mean(retest_thresholds)];
N = length(sham_thresholds);
err1 = std(sham_thresholds)/sqrt(N);
err2 = std(retest_thresholds)/sqrt(N);
CI95 = tinv([0.025 0.975], N-1);
err1CI95 = bsxfun(@times, err1, CI95(:));
err2CI95 = bsxfun(@times, err2, CI95(:));
err = [err1CI95(2) err2CI95(2)];

figure();
errorbar(x,y,err);
xlim([-0.25 1.25]);
ylim([-1 10]);
xlabel('SR condition')
ylabel('Auditory Threshold (dB)');
title('Plot to show re-test values');
hold on;

for i = 1:length(sham_thresholds)
    plot(x, [sham_thresholds(i) retest_thresholds(i)], 'color', [0.8 0.8 0.8]);
end

improvement = retest_thresholds - sham_thresholds;
[R,P] = corrcoef(sham_thresholds,-1*improvement)

figure()
D = sham_thresholds;
F = improvement;  
scatter(D,F);
p = polyfit(D,F,1);
f = polyval(p,D);
hold on
plot(D,f,'-r');
xlabel('Re-Test Sham Auditory Threshold (dB)');
ylabel('Absolute Improvement in Auditory Threshold (dB)');
titlename = ['Improvements against Re-Test Sham, corrP = ' num2str(P(1,2))];
title(titlename);
