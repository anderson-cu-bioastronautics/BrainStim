% Bubble plot and staircase function

function [x,prob_bi,prob_hat] = plotSBC(stim,cor,lambda_hat,mu_bi,sigma_bi,mu_hat,sigma_hat,guess_rate,titlestr,xstr,ystr,clusterSize)

% x and y to plot the cdf
x = linspace(min(stim)-0.1,max(stim)+0.5,400);

% Cumulative normal distribution of the fitted psychometric curve
prob_bi = guess_rate + (1-guess_rate-0)*cdf('norm', x, mu_bi, sigma_bi);

% Cumulative normal distribution of the fitted psychometric curve with
% a guess rate differing from 0.5 and non-zero lapse rate
prob_hat = guess_rate + (1-0.5-lambda_hat)*cdf('norm', x, mu_hat, sigma_hat);

figure;
%% The staircase (S)
subplot(1,3,1); hold on;
for i = 1:length(stim)
    if cor(i)
        plot(i,stim(i),'ro');
    else
        plot(i,stim(i), 'rx');
    end
end
hold off;

%     plot((1:length(stim)),stim,'r.','MarkerSize',20);
title(titlestr);
xlabel('Presentation');
ylabel('Stimulus');
grid minor;

%% The bubble (B)
subplot(1,3,2);
plot(x,prob_bi,'g','LineWidth',1.5);
hold on;
plot(x,prob_hat,'b--','LineWidth',1.25);
% Average the Y outcomes at each unique X value
[Xunique, ix, ixu] = unique(stim);
Punique = zeros(1,length(Xunique));
Lunique = zeros(1,length(Xunique));
for k = 1:length(Xunique)
    YatXunique = cor(ixu == k); % find the Y outcomes for the jth unique X value
    Lunique(k) = length(YatXunique);    % find the number of trials for the ith unique X value
    Punique(k) = mean(YatXunique);  % find the probability at the ith unique X value
end
for k = 1:length(Xunique)
    % find the marker size with a max of 12 and a min of 4
    if max(Lunique) == min(Lunique)
        msize = 8;
        
    else
        msize = (12-4)/(max(Lunique)-min(Lunique)) * (Lunique(k)-min(Lunique)) + 4;
    end
    plot(Xunique(k), Punique(k), 'ko', 'MarkerSize', msize);
end
title(strcat(titlestr,' with Bubble Plot'));
xlabel(xstr);
ylabel(ystr);
grid on;
axis([min(stim)-0.2 max(stim)+0.2 0 1]);
legend('Binary Psychometric Curve','Lapse Psychometric Curve','Relative Data Points','Location','southeast');

%% The cluster (C)

% How many clusters are we plotting?
numClusters = ceil(length(stim)/clusterSize);

% Variables to hold the values to plot
X = zeros(1,numClusters);
Y = zeros(1,numClusters);

[sortStim, isort] = sort(stim);
sortCor = cor(isort);
% For loop for all the clusters
for i = 0:numClusters-1
    X(i+1) = mean(sortStim(i*clusterSize+1:(i+1)*clusterSize));
    Y(i+1) = mean(sortCor(i*clusterSize+1:(i+1)*clusterSize));
end
if mod(length(stim),clusterSize) ~= 0
    X(numClusters) = mean(sortStim((i+1)*clusterSize):length(stim));
    Y(numClusters) = mean(sortCor((i+1)*clusterSize):length(stim));
end

% This Y should NOT just be the psychometric curve!
% It should be the average of the values for each cluster
% Y = guessRate + (1 + guessRate - 1 - lambda)*cdf('normal',X,mu,sigma);

subplot(1,3,3);
hold on;
plot(x,prob_bi,'g','LineWidth',1.5);
plot(x,prob_hat,'b--','LineWidth',1.25);
plot(X,Y,'b.','MarkerSize',30);
title(strcat(titlestr,' with Cluster Plot'));
xlabel(xstr);
ylabel(ystr);
grid on;
axis([min(stim)-0.2 max(stim)+0.2 0 1]);
legend('Binary Psychometric Curve','Lapse Psychometric Curve','Relative Data Points','Location','southeast');
end