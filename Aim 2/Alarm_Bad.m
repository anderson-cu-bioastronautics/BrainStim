function [] = Alarm_Bad(start_time)
x = linspace(0,15+start_time,(15+start_time)*11000);%[0:1/11000:1];
y =-3*ones(1,start_time*11000); 
y=[y linspace(-2,-1.75,15*11000)];%;[-3:1/11000/3:0];
WarnWave = [sin(2*pi*1000*x)];
WarnWave = WarnWave+y;
%plot(WarnWave)
%audiowrite('Alarm_10s_Constant.wav',WarnWave,10800);
sound(WarnWave,10800)

