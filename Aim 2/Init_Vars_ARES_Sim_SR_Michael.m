%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defines the initial conditions and constant variables for
% the flight simulator and has been adapted for the LM USER Experiment.
% Elliott Davis 
% Edited: 4/30/18
% Carlos Pinedo
% Changed to accomodate time delay study.
% Edited: 4/22/18
% Nicholas Miller
% Added failure functionality for LM USER Experiment.
% Edited: 
% Daniel Gutierrez Mendoza
% Randomized and reset initial coditions in each scenario.
% safetAdded 'PreInitVars' function, etc to accomodate double-blind.
% Edited: 6/7/2019 
% Daniel Gutierrez Mendoza
% Menu option Training/Testing.
% Edited: 6/12/2019
% Joshua Seedorf
% Complete rework for realtime ALD with new maps
% 5/5/2020
% Edited: 09/05/2020
% Generalized for ARES Sim use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
global R0val 
rng('shuffle');
% Menu for deciding whether we are testing or training. (Is this necessary
% since they are only going to be training a few times and then every trial
% afterwards will be testing.

% Check that Post Hazards Ran
if exist('safety','var')
    if ~safety
        error('Run Post Hazards or clear variables for a new Subject!')
    end
end

clearvars -except Subject safety trials_done use_ald modes_available...
          fuel_order trial_order Avail_Mode
clc

%% Select Trial 
%trial_num=input('Enter Trial Number: \n');
% if ~exist('trials_done','var') %If first initialization
%     safety = true;
%     
%     trials_done = [];
%     
%     %TEMPORARY SUBJECT FOR TESTING
%     Subject = input('ENTER SUBJECT NUMBER: ');
% 
%     %TEMPORARY MODE AVAILABLE
%     Avail_Mode = 0; %[0 = RCAH , 1 = MIPC , 2 = AIPC]
%     
%     % Randomize the Trial order
%     trial_order_1 = randperm(8);
%     
%     trial_order_2 = randperm(8);
%     
%     % Prevent back to back trials
%     run = false;% COMMENTED WHLE WORKING WITH MAPS
%     for tt = 1:3
%         run = run + ~isempty(find(trial_order_1(5+tt) == trial_order_2(1:3), 1));
%     end
%     
%     while run
%         trial_order_2 = randperm(8);
%         run = false;
%         for tt = 1:3
%             run = run + ~isempty(find(trial_order_1(5+tt) == trial_order_2(1:3), 1));
%         end
%     end
%     
%     trial_order = [trial_order_1,trial_order_2] + 12*ones(1,16);
%     trial_order = [1:12,trial_order];
%     
%    clear trial_order_1 trial_order_2
%     
%  
%     
% 
% end
trials_done=[];
% trial_order=[1:12,1:20]; % Remove after MAPS finshed 
trial_order=[1:9,1:36]; % Remove after MAPS finshed 

% for x = 1:12
for x = 1:9
    if isempty(trials_done(trials_done == x))
        if x < 10
            str = sprintf('TRAINING 0%d [ ]',x);
        else
            str = sprintf('TRAINING %d [ ]',x);
        end
    else
        if x < 10
            str = sprintf('TRAINING 0%d [X]',x);
        else
            str = sprintf('TRAINING %d [X]',x);
        end
    end
    list(x) = {str};
end
% for n = 1:32

for n = 1:36
    if isempty(trials_done(trials_done == n+x))
        if n < 10
            str = sprintf('Trial 0%d [ ]',n);
        else
            str = sprintf('Trial %d [ ]',n);
        end
    else
        if n < 10
            str = sprintf('Trial 0%d [X]',n);
        else
            str = sprintf('Trial %d [X]',n);
        end
    end
    list(n+x) = {str};
end

[trial_num,tf] = listdlg('ListString',list);

trials_done(end+1) = trial_num;

%Initial mass of fuel (either 50 or 70 slugs)
m_fuel0 = 70; 


scenario_num = trial_order(trial_num);


if trial_num <= 9
    fprintf('Running TRAINING %d...\n', trial_num)
    train = true;
  
else
    fprintf('Running Trial %d...\n', trial_num - x)
    train = false;
end

%% Set constants
% Vehicle/Engine Characteristics
Isp = 311;      % sec, specific impulse of descent engine thruster fuel
g0 = 9.81 * 3.28084;   % ft/s^2, gravitational constant used for Isp calc
m_Lander = 493; % slugs

% Initial Position (center of map)
x0 = 12*2000/2; %ft
y0 = 2*12*2000/5; %ft

% Planetary Characteristics
gBody = 1.622 * 3.28084;   % ft/s^2, lunar gravity

% Control Constants
tauV = 8;       % seconds, time constant of velocity guidance
tauh = 25;      % seconds, time constant of height guidance
tauT = 1.5;     % seconds, time constant of thrust control

% Define the low gate and terminal descent values
hdotTD = -3;    % ft/s
hLG = 500;      % ft
hdotLG = -16;   % ft/s
hTD = 150;      % ft

a = (hdotLG - hdotTD)/(hLG-hTD);
b = -a*hTD + hdotTD;
C = -b/a;
k = hLG-C;
q = (2*a)/log(k/(hTD-C));



%% Failure Variables, left from USER
fail_type=0; % 1=Pitch; 2=roll;3=cont rev;4=stuck thruster;5=fuel gauge 6=roll triag 
fail_time=1300;
fail_detect=0;
fail_time_s = fail_time/30; %fail time in seconds for gradual failures


quad_s = 1;
    
%% For Aim 2
if trial_num<10
    nhname = sprintf('SR_MAPS\\maps\\map_0%d_nh.png',trial_num); %Image no hazards
    yhname = sprintf('SR_MAPS\\maps\\map_0%d_yh.png',trial_num); %Image with hazards (hidded for training)
    dtname = sprintf('SR_MAPS\\maps\\map_0%d_data.mat',trial_num); %Map Data
    hsname = sprintf('SR_MAPS\\maps\\map_0%d_haz.png',trial_num); %HSI map
else 
    nhname = sprintf('SR_MAPS\\maps\\map_%d_nh.png',trial_num); %Image no hazards
yhname = sprintf('SR_MAPS\\maps\\map_%d_yh.png',trial_num); %Image with hazards (hidded for training)
dtname = sprintf('SR_MAPS\\maps\\map_%d_data.mat',trial_num); %Map Data
hsname = sprintf('SR_MAPS\\maps\\map_%d_haz.png',trial_num); %HSI map
    
end   


image_nh = imread(nhname);image_yh = imread(yhname);image_hs = imread(hsname);
[size_mat1,size_mat2,size_mat3]=size(image_nh);
if size_mat1~=1124 || size_mat2~=1124 
    image_nh=image_nh(9:1132,1:1124,1:3);
    image_yh=image_yh(9:1132,1:1124,1:3);
    image_hs=image_hs(9:1132,1:1124,1:3);
end 
load(dtname)

mapping = 24;
ftperpix = mapping*2000/length(image_nh);

% After conditional write to current chosen landing points
csvwrite('PntsMatFt.csv',LPs)
csvwrite('PntsMatPix.csv',round(LPs/ftperpix))
csvwrite('SciPntsMatFt.csv',SPs)
csvwrite('SciPntsMatPix.csv',round(SPs/ftperpix))
imwrite(image_nh,'ALD_Map_No_Hazzards.png')
imwrite(image_yh,'ALD_Map_Hazzards.png')
imwrite(image_hs,'VRWorlds\\ALD_Map_HSI.png')

% Bilimoria Initial Landing Point
xLP0 = x0; %Initialize starting point of the lander
yLP0 = y0 + 1350; 

%% Vehicle Constants
h0 = 500; % Initial height
R0 = sqrt((y0-yLP0)^2 + (x0-xLP0)^2); % ft, initial range distance
R0Val = R0; %Save R0 as separate value
R_R0_tol = exp(1)^2 - 1.01; %Tolerance between R and initial R
hLP = 0 ; %Height of the landing point (For these simulations the terrain height always=0)
hTD = hLP + 150; % Height at Terminal Descent
m0 = m_Lander + m_fuel0; % Slugs, initial wet mass of the vehicle
VN0 = 60 ; %Initial Northern velocity of spacecraft
VE0 = 0 ; %Initial Eastern velocity of spacecraft
theta0 = 16 ; %Initial pitch angle
phi0 = 0 ; %Initial roll angle
hdot0 = -10 ; %Initial rate of descent
%fprintf('WARNING: \ng = %d, pitch = %d, and roll = %d for troubleshooting PFD\n',gBody,theta0,phi0)

%% Define the tolerances used throughout the program
fzero_search_tolerance = 0.05;      % percentage of ft used. implemented via toler = 2.0*tol*max(abs(b),1.0);
terrain_step_for_zero_check = 10;   % ft - I dont think this being implemented anywhere
MP_step = 10;   % ft
odeRelTol = 0.01;   % Default: 10^-3
odeAbsTol = 10^-3;  % Default: 10^-6
% tolerance for what is an acceptable distanceif h0 < 0
landing_range_error = 50;   % ft
% pts used to check if the zero fnd is the true first zero
ALD1_backcheck_pts = 5;

tolerance_vec = [fzero_search_tolerance,terrain_step_for_zero_check, MP_step,landing_range_error,ALD1_backcheck_pts,odeRelTol,odeAbsTol];

%% Define tolerances used in guidance and trajectory calculations
% tolerance for how far down range is acceptable to be in a hover until TD
Rlim = 50;  % ft
% tolerance for what is an acceptable value for R/R0 < e^2 - tol
% R_R0_tol = exp(1)^2 -1; % this value was selected abritrily 
% my thought process was that the orginial Bilimoria set up wanted pilots
% to fly towards the landing site so there was no need to move backwards.
% Subtracting 1 means that there is a tighter constraint on R being greater
% than R0 therefore it is more likely to generate a new trajectory as
% opposed to subtracting say 6. (As R_R0_tol decreases, the constraint on
% the range R exceeding the initial range R0 becomes more relaxed and will
% likely only generate a new guidance trajectory once or twice in a sim)
max_thrust = 10000; %lbf, defined by Bilimoria
thrust_u = .6;  % 60%, defined by Bilimoria
thrust_l = .1;  % 10%, defined by Bilimoria
max_phi = 45;   % deg, defined by Bilimoria
max_theta = 45; % deg, defined by Bilimoria
sim_tol = [Rlim,R_R0_tol,max_thrust,thrust_u,thrust_l,max_phi,max_theta];

%% Run Trajectory Simulation
% Define the initial height of the vehicle relative to the height of the LS
h0_sim = h0 - hLP;
% define initial state
X0 = [y0, x0, VN0, VE0, theta0, phi0, h0_sim, hdot0, m0];
% define the gain for the model pilot
K = 1;
vars = [R0 hTD tauV tauh gBody tauT Isp g0 a b C k q m_Lander K mapping];
R0val = vars(1);
% IPC interval
interval = 10;
mp_guess = 0;


%% Tactile Stuff
set_Ex(); % sets  global variables in extrinsic matlab 
buzz_time_1=3; 
alarm_time_1=randi([0 10]);
buzz_time_2=15;
%% Constants for simulink
Time1 = 0;
Time2 = 0;
Time3 = 0;


input_vars = [0,R0, gBody, a, b, C, k, q, xLP0, yLP0, hTD, hLP];
constant_vars = [0,tauV, tauT, tauh, g0, Isp, m_Lander, R_R0_tol,interval];
Initial_conditions = [Time3, m_fuel0, x0, y0, VE0, VN0,...
                      theta0, phi0, h0, hdot0, mp_guess, fail_type, ...
                      fail_time,fail_time_s,fail_detect,quad_s,buzz_time_1,alarm_time_1,buzz_time_2];  
                  
                  
Avail_Mode=0;                  
%% Tactile Stuff
set_Ex(); % sets  global variables in extrinsic matlab 
ping=menu('Do you require to ping tactor?', 'Yes','No');

if ping==1
    CallTactile(.7,400,1);pause(1);CallTactile(.7,400,1); % pings tactle device
end
safety = false;
                  

% %% Calculates best landing point 
% % Calculate the mean science point
% Xffeet=Xf; %[ft]
% Yffeet=Yf; %[ft]
% numPoints = 2000;
% feetPerPoint = 24000/numPoints;
% Xf = Xf/feetPerPoint; %[points]
% Xf = round(Xf);
% Yf = 2000-Yf/feetPerPoint;
% Yf = round(Yf); %[points]
% SP_mean = [24000 - mean(SPs(2,:)) ; 24000 - mean(SPs(1,:))];
% 
% % Calculate the distance for each landing point
% for n = 1:6
%     LP = [LPs(1,n);24000 - LPs(2,n)]; %[feet]
%     dist_SP_NH(n) = norm(SP_mean - LP); %[feet]
%     dist_LP(n) = norm(LP - [Xffeet;24000-Yffeet]); %[feet]
%     
%     % Eliminate Hazardous Points
%     LP = LP/feetPerPoint; %[points]
%     X_LP = (LP(1) - round(lander_size/(2*feetPerPoint)):(LP(1) + round(lander_size/(2*feetPerPoint)) - 1)); %[points]
%     Y_LP = (2000 - LP(2) - round(lander_size/(2*feetPerPoint))):(2000 - LP(2) + round(lander_size/(2*feetPerPoint)) - 1); %[points]
%     
%     LP_hazards = 100*sum(sum(hazards(Y_LP,round(X_LP))))/(lander_size/feetPerPoint)^2; %[percent]
% 
%     if LP_hazards > 25
%         dist_SP(n) = Inf;
%     else
%         dist_SP(n) = dist_SP_NH(n);
%     end
% end
% % Find the best
% best_LP = find(dist_SP == min(dist_SP));
                  
 