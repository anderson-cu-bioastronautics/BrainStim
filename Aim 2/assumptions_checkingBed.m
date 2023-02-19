function stu_res = assumptions_checking(p,tbl,stats,terms,Performance,Condition,Subjects)
% clear all; close all; clc;
% 
% %[p,tbl,stats,terms] = anovan(hs_pns(:,3),{scene,time_short,subj},'model',[1,0,0; 0,1,0; 0,0,1; 1,1,0],'random',3,'varnames',{'Scene','Time','Subj'});
% Data = readtable('LOTNoiseAccuracyData.csv');
% Subjects = Data(:,1);
% Meas = table([1 2 3 4]','VariableNames',{'Conditions'});
% rm = fitrm(Data,'Sham-MMSR~Subjects','WithinDesign',Meas);
% ranova(rm)
%% Diagnostics

res = stats.resid; %Response_all(1:36,3);
StatTable = [Condition Performance];
% 
% res = residuals(dat_compTime_stats);
% y_pred = mean(data);
% 
% res = zeros(length(data),1);
% for i = 1:length(data)
%     res(i) = data(i) - y_pred;
% end
x1 = [1:length(res)]';

% check for homoscedasticity
figure()
subplot(2,1,1)
scatter(x1,res) 
title('Residuals')
ylabel('Residuals')
hold on
plot(x1,zeros(length(x1),1),'r')


% studentized residual plot
%MSE = tbl{5,5};
MSE = tbl{4,5};
stu_res = res./sqrt(MSE);
subplot(2,1,2)
scatter(x1,stu_res)
title('Semi-Studentized Residuals')
ylabel('Semi-Studentized Residuals')
hold on
plot(x1,zeros(length(x1),1),'r')

% % check for independence
% % by task
% for i = 1:length(Task_num)
%     if Task_num(i) == 1
%         s1(i) = res(i);
%     elseif Task_num(i) == 2
%         s2(i) = res(i);
%     elseif Task_num(i) == 3
%         s3(i) = res(i);
%     elseif Task_num(i) == 4
%         s4(i) = res(i);
%     elseif Task_num(i) == 5
%         s5(i) = res(i);
%     elseif Task_num(i) == 6
%         s6(i) = res(i); 
%     elseif Task_num(i) == 7
%         s7(i) = res(i);        
%     end
% end
% figure()
% subplot(1,2,1)
% bar([1;2;3;4;5;6;7],[mean(abs(s1));mean(abs(s2));mean(abs(s3));mean(abs(s4));mean(abs(s5));mean(abs(s6));mean(abs(s7))])
% title('Avg Residuals by Task')
% xlabel('Task')
% ylabel('Avg Absolute Residual')
% xticks([1:7])
% xvalues = {'DSST';'LOT';'MPT';'MRT';'NBACK';'PVT';'VOLT'}
%     xticklabels(xvalues);
%     xlim([0 8])
% subplot(1,2,2)
% scatter(Task_num,res)
% title('Residuals by Task')
% xlabel('Task')
% ylabel('Residual')
% xticks([1:7])
% xvalues = {'DSST';'LOT';'MPT';'MRT';'NBACK';'PVT';'VOLT'}
%     xticklabels(xvalues);
% xlim([0 8])
% hold on
% plot(linspace(0,8,length(x1)),zeros(length(x1),1),'r')

% by condition
for i = 1:length(Condition)
    if Condition(i) == 1
        t1(i) = res(i);
    elseif Condition(i) == 2
        t2(i) = res(i);
    elseif Condition(i) == 3
        t3(i) = res(i);
    elseif Condition(i) == 4
        t4(i) = res(i);        
    end
end
figure()
subplot(1,2,1)
bar([1;2;3;4],[mean(abs(t1));mean(abs(t2));mean(abs(t3));mean(abs(t4))])
title('Avg Residuals by Condition')
xlabel('Condition')
ylabel('Avg Absolute Residual')
subplot(1,2,2)
scatter(Condition,res)
xlim([0 5])
title('Residuals by Condition')
xlabel('Condition')
ylabel('Residual')
hold on
plot(linspace(0,5,length(x1)),zeros(length(x1),1),'r')

% % by subject
% for i = 1:length(subj)
%     if subj(i) == 1
%         sub1(i) = res(i);
%     elseif subj(i) == 2
%         sub2(i) = res(i);
%     elseif subj(i) == 3
%         sub3(i) = res(i);
%     elseif subj(i) == 4
%         sub4(i) = res(i);
%     elseif subj(i) == 5
%         sub5(i) = res(i);
%     elseif subj(i) == 6
%         sub6(i) = res(i);
%     end
% end
% figure()
% subplot(1,2,1)
% bar([1;2;3;4;5;6],[mean(abs(sub1));mean(abs(sub2));mean(abs(sub3));mean(abs(sub4));mean(abs(sub5));mean(abs(sub6))])
% title('Avg Residuals by Subject')
% xlabel('Subject')
% ylabel('Avg Absolute Residual')
% subplot(1,2,2)
% scatter(subj,res)
% xlim([0 7])
% title('Residuals by Subject')
% xlabel('Subject')
% ylabel('Residual')
% hold on
% plot(linspace(0,7,length(x1)),zeros(length(x1),1),'r')

% check for normal distribution
figure()
subplot(2,1,1)
histogram(res) 
title('Histogram of Residuals')
subplot(2,1,2)
histogram(stu_res)
title('Histogram of Semi-Studentized Residuals')
xlabel('Residuals')

% normal probability plot
figure()
normplot(res)
title('Normal Probability Plot of Residuals')
ylabel('Residual')

[h,p] = adtest(stu_res)
vartestn(Performance,Condition)

god_hates_me = 1;

% % 
% figure()
% plotResiduals()