;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Hazards Realtime
%   Saves data for the realtime ALD experiement in the AReS.
%   Joshua Seedorf 5/14/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
Subject = input("Subject Number : ");
%% Save the raw workspace
save('rawtmp.mat')

%% Create Directories
if ~exist('Subject_Data','dir') % "Subject_Data" Folder
    mkdir Subject_Data
end

%Subject
if Subject >9
    Subject_str = num2str(Subject);
else
    Subject_str = ['0',num2str(Subject)];
end 
subFolder = ['Subject_Data//Subject_',Subject_str];

if ~exist(subFolder,'dir') % "Subject_##" Folder
    mkdir(subFolder);
end

%% Keep important variables
%   THIS COMMAND IS REPEATED AT THE END OF THIS SCRIPT, COPY ANY CHANGES
clearvars -except Subject trial_num scenario_num use_ald hazards abort_press... 
    fail_detection fuel h hdot pitch roll Time VE VN xpos ypos valid_abort...
    train trials_done xLP yLP Avail_Mode thetaG phiG fuel_order trial_order...
    subFolder Subject_str safety LPs SPs RealTime LP_button IPC_button alarm_time_1 Auditory_TIME_BUTTON...
    m_Lander
     
%% Get Vehicle Information
Xf=xpos.data(end);
Yf=ypos.data(end);
Xf_f=Xf;
Yf_f=Yf;

abort = valid_abort.data(end);

if ~train
    trial_num = trial_num - 12;
end

if trial_num > 9
    trial_num_str = num2str(trial_num);
else 
    trial_num_str = ['0',num2str(trial_num)];
end

if scenario_num > 9
    scenario_num = num2str(scenario_num);
else 
    scenario_num = ['0',num2str(scenario_num)];
end



Xffeet=Xf; %[ft]
Yffeet=Yf; %[ft]
numPoints = 2000;
feetPerPoint = 24000/numPoints;
Xf = Xf/feetPerPoint; %[points]
Xf = round(Xf);
Yf = 2000-Yf/feetPerPoint;
Yf = round(Yf); %[points]

%% Decide if the vehicle landed in a hazard or not. 
lander_size = 50; %[ft]

% Pull the nearest 7x7 area to evaluate
X_v = (Xf - 3):(Xf + 3); %[points]
Y_v = (2000 - Yf - 3):(2000 - Yf + 3); %[points]

% Evaluate the lander coverage of every point
weighted_hazards = zeros(6);

for nx = 1:7
    for ny = 1:7
        % 0 = none, 1 = hazards
        hazard_level = hazards(Y_v(ny),X_v(nx));
        
        % Determine the physical area of this point which is covered by the lander
        dx = lander_size/2 - (abs(Xffeet - X_v(nx)*feetPerPoint) - feetPerPoint/2); %[feet]
        dy = lander_size/2 - (abs(Yffeet - Y_v(ny)*feetPerPoint) - feetPerPoint/2); %[feet]
        
        % Limit the maximum dimension to the size of the point and avoid
        % negatives
        if dx > feetPerPoint
            dx = feetPerPoint;
        elseif dx < 0
            dx = 0;
        end
        if dy > feetPerPoint
            dy = feetPerPoint;
        elseif dy < 0 
            dy = 0;
        end
        
        % Calculated the weighted hazards at this point
        weighted_hazards(ny,nx) = ((dx*dy)/lander_size^2) *hazard_level;
    end
end

% Sum the weighted hazards to get the percent_hazards
if ~train
    percent_hazards = sum(sum(weighted_hazards))*100;
elseif trial_num >= 5
    percent_hazards = sum(sum(weighted_hazards))*100;
else
    percent_hazards = 0;
end


%% Roll and Pitch angles
Roll_angle = roll.data(end);
Pitch_angle = pitch.data(end);
Tilt = norm([Roll_angle,Pitch_angle]); % Max 6 deg

%%  Velocity Calculation
Speed = norm([VN.data(end),VE.data(end)]); % Max 3 ft/s

%% Descent Rate and Fuel
Des_Vel = abs(hdot.data(end)); % Max 4 ft/s
Fuel    = fuel.data(end);
height  = h.data(end);

fuel_crash = Fuel <= 0 && height >= 3; % Check if fuel ran out more than three feet above the ground

%% Final Crash Determination
clc
if ~abort
    if percent_hazards > 25 || Tilt > 6 || Speed > 4 || Des_Vel > 8 || fuel_crash
        crash = true;
        fprintf('CRASHED!\n')
        if percent_hazards > 25
            fprintf('Landed %d%% in hazards\n',round(percent_hazards))
        end
        if Tilt > 6
            fprintf('Tilt of %d degrees exceeded the 6 degree limit\n',ceil(Tilt))
        end
        if Speed > 4
            fprintf('Speed of %d ft/s exceeded the 4 ft/s limit\n',ceil(Speed))
        end
        if Des_Vel > 8
            fprintf('Descent velocity of %d ft/s exceeded the 8 ft/s limit\n',ceil(Des_Vel))
        end
        if fuel_crash
            fprintf('Ran out of fuel at %d ft\n',ceil(height))
        end
    else
        crash = false;
        fprintf('LANDED!\n')
        if trial_num >= 5 && train
        fprintf('Landed %d%% in hazards\n',round(percent_hazards))
        end
    end
else
    crash = false;
    fprintf('ABORTED!\n')
end




%% Calculates best landing point 
% Calculate the mean science point
SP_mean = [24000 - mean(SPs(2,:)) ; 24000 - mean(SPs(1,:))];

% Calculate the distance for each landing point
for n = 1:6
    LP = [LPs(1,n);24000 - LPs(2,n)]; %[feet]
    dist_SP_NH(n) = norm(SP_mean - LP); %[feet]
    dist_LP(n) = norm(LP - [Xffeet;24000-Yffeet]); %[feet]
    
    % Eliminate Hazardous Points
    LP = LP/feetPerPoint; %[points]
    X_LP = (LP(1) - round(lander_size/(2*feetPerPoint)):(LP(1) + round(lander_size/(2*feetPerPoint)) - 1)); %[points]
    Y_LP = (2000 - LP(2) - round(lander_size/(2*feetPerPoint))):(2000 - LP(2) + round(lander_size/(2*feetPerPoint)) - 1); %[points]
    
    LP_hazards = 100*sum(sum(hazards(Y_LP,round(X_LP))))/(lander_size/feetPerPoint)^2; %[percent]

    if LP_hazards > 25
        dist_SP(n) = Inf;
    else
        dist_SP(n) = dist_SP_NH(n);
    end
end
% Find the best
best_LP = find(dist_SP == min(dist_SP));

% Find the closest landing point
actual_LP = find(dist_LP == min(dist_LP));


distance2Best_LP = dist_LP(best_LP);
distance2MeanSP = norm(SP_mean - [Xffeet;24000-Yffeet]);
differenceFromOptimal = distance2MeanSP - dist_SP(best_LP);

% Check Optimization
%   10% Buffer on "optimal" range
if distance2MeanSP < 450 + dist_SP_NH(best_LP)
    % Landed in an optimal zone
    optimal = 1;
else
    % Landed in a sub-obtimal zone
    optimal = 0;    
end

%% Money conditional 
if ~train
    fprintf('\n')
    if crash && ~abort
        bucks = -5;
        fprintf('You CRASHED   (And lost $5)! \n')
    elseif ~crash && ~abort     
        if ~optimal
            fprintf('You did not optimize your landing zone!\n   You earned $1\n')
            bucks = 1;
        elseif optimal
            fprintf('You optimized your landing zone!\n   You earned $2\n')
            bucks = 2;
        end
    end 
    if abort 
        bucks = 0;
        abort=true;
        fprintf('You ABORTED, you get zero dollars... better than -5\n')
        crash=false;
    end 
else
    if ~optimal && ~crash && ~abort
        fprintf('You did not optimize your landing zone!\n')
        
    elseif optimal && ~crash && ~abort
        fprintf('You optimized your landing zone!\n')
        
    end
    abort=0;
end 

if train==true
    bucks=0;
    crash = false;
end 

%% Calculates the total bonus up to now
if ~train && ~safety
    if exist('payment.mat','file')
        load('payment.mat');
        if length(total_payment) >= Subject
            total_payment(Subject) = total_payment(Subject) + bucks;
        else
            total_payment(Subject) = bucks;
        end
    else
        total_payment(Subject) = bucks;
    end

    if total_payment(Subject) < 0 && trial_num == 16
        total_payment(Subject) = 0;
    end
    
    fprintf('\nYour bonus is now $%d\n',total_payment(Subject))
    save('payments.mat','total_payment')
end 

%% Landing Point Selection
LP_selection(1) = 0; %Bilmoria Point
for step = 2:length(LP_button.data)
    if LP_button.data(step) ~= 0
        LP_selection(step) = LP_button.data(step);
    else
        LP_selection(step) = LP_selection(step-1);
    end
end
LP_selection = timeseries(LP_selection,LP_button.time);

%% Number of Redesignations
% Include first from Bilmoria
num_selections = 0;
for step = 2:length(LP_selection.data)
    if LP_selection.data(step) ~= LP_selection.data(step-1)
        num_selections = num_selections +1;
    end
end

%% Reaction Time
% Test times
times = [14,26,39,57,73,93,106,126,143,155,171,189,205,223,240]; %[seconds]
% Max duration for each test 
dt = 10; %[seconds]

reaction_times = zeros(1,length(times));

for event = 1:length(times)
    % Events are only valid if the entire 10 second duration occurs during
    % the simulation
    if (times(event) + dt) < Time.data(end)
        % Find the press
        tmin = find(fail_detection.time >= times(event),1);
        tmax = find(fail_detection.time >= times(event)+dt,1);
        tspan = [tmin:tmax];
        detect = find(fail_detection.data(tspan) == 1,1);
        
        % Calculate Reaction Time
        if ~isempty(detect)
            reaction_times(event) = fail_detection.time(tspan(detect)) - times(event);
        else
            % Missed the event
            reaction_times(event) = Inf; % should this be 10 sec? 
        end
    else
        % Did not occur
        reaction_times(event) = -1;
    end
end
hits = find( reaction_times(reaction_times ~= -1) ~= Inf);
misses =  count(num2str(reaction_times),'Inf');
mean_rt = mean(reaction_times(hits));
%   fprintf('\nReaction Time Average    %.3f sec\n',mean_rt);
%   fprintf('Number of Missed         %d\n', misses);
%% Guidance Error
% Taken for pitch and roll command and a third for combined angle
pitch_error = zeros(1,length(thetaG.data));
roll_error = zeros(1,length(phiG.data));
for step = 1:length(thetaG.data)
    pch = pitch.data(find(pitch.time >= thetaG.time(step),1));
    rol = roll.data(find(roll.time >= phiG.time(step),1));
    
    if isempty(pch) || isempty(rol)
        pch = pitch.data(end);
        rol = roll.data(end);
    end
    
    pitch_error(step) = pch - thetaG.data(step);
    roll_error(step) = rol - phiG.data(step);
end


    total_error = sqrt(pitch_error.^2 + roll_error.^2);

    rms_pitch = rms(pitch_error);
    rms_roll  = rms(roll_error);
    rms_total = rms(total_error);
    FPM = 100/(1+rms_total);
    %fprintf('Flying Performance Rating    %.1f out of 10\n',FPM);
    
%% Crash - abort - land
    
%tspan = [Time.Data(end) (Time.Data(end)+1000)]; 
tspan = [0 1000];

Initial = [xpos.Data(end), ypos.Data(end), VN.Data(end), VE.Data(end), pitch.Data(end), roll.Data(end), h.Data(end), hdot.Data(end), (m_Lander + fuel.Data(end))];

% ODE45 
opts = odeset('Events', @stopping_point);
[tOUT,xOUT] = ode45(@trajectory_sim, [tspan],Initial, opts);


xPos_ode = xOUT(:,1);% x position output from ode
yPos_ode = xOUT(:,2); % y position output from ode
VN_ode = xOUT(:,3); 
VE_ode = xOUT(:,4);
pitch_ode = xOUT(:,5);
roll_ode = xOUT(:,6);
height_ode = xOUT(:,7);
hDot_ode = xOUT(:,8);
mass_ode = xOUT(:,9);
time_ode = tOUT(:,1);

% Finding the distance they were projected to be able to travel
initialPos = sqrt(xPos_ode(1)^2 + yPos_ode(1)^2);
endPos = sqrt(xPos_ode(end)^2 + yPos_ode(end)^2);
distProjected = abs(endPos - initialPos);


if abort == 0 % If they did not abort
    if crash == 0 % if they did not crash
        CAL_Score = 1; % then they landed.
    elseif crash == 1 % Otherwise, they did crash.
        CAL_Score = 0;
    end
    
elseif abort == 1 % if they aborted
    % Find the distance to any landing point from the ODE output 
    for i = 1:length(LPs(1,:))
        for j = 1:length(xPos_ode)
            xDistToLP(j,i) = LPs(1,i) - xPos_ode(j);
            yDistToLP(j,i) = LPs(2,i) - yPos_ode(j);
            distToLP(j,i) = sqrt(xDistToLP(j,i)^2 + yDistToLP(j,i)^2); % distance to every landing point when they run out of fuel
        end
    end
    
    % indexing the minimum distance to every landing point
    for i = 1:length(distToLP(1,:))
        [valMinDist_ode(i),idxMinDist_ode(i)] = min(distToLP(:,i));
    end
    
    % indexing when that happens and to which landing point
    [minIdxODE,idxMinIdxODE] = min(idxMinDist_ode);

    if valMinDist_ode(idxMinIdxODE) <= distProjected % if projected distance is greater than min dist from LP, they could have landed
        CAL_Score = 0.333;
    elseif valMinDist_ode(idxMinIdxODE) > distProjected % if projected distance is less than min dist from LP, they couldn't have landed
        CAL_Score = 0.667;
    end

    
end
%% Save Data
savedata = menu('Save Data?','YES','No');
if savedata == 1
    
    
    if train
        filename = [subFolder,'//','Subject_',Subject_str,'_Training_',trial_num_str,'.mat'];
        rawfilename = [subFolder,'//','Subject_',Subject_str,'_Training_',trial_num_str,'_RAW.mat'];

        copy_num = 0;
        if exist(filename,'file')
            rep = menu('That File Alread Exists!','Replace the file','Save as a new file');
            if rep == 2
                while exist(filename,'file')
                    copy_num = copy_num + 1;
                    copy_txt = sprintf('(%d)',copy_num);
                    filename = [subFolder,'//','Subject_',Subject_str,'_Training_',trial_num_str,'_',copy_txt,'.mat'];
                    rawfilename = [subFolder,'//','Subject_',Subject_str,'_Training_',trial_num_str,'_',copy_txt,'_RAW.mat'];
                end
            end
        end
    else
        filename = [subFolder,'//','Subject_',Subject_str,'_T',trial_num_str,'.mat'];
        rawfilename = [subFolder,'//','Subject_',Subject_str,'_T',trial_num_str,'_RAW.mat'];

        copy_num = 0;
        if exist(filename,'file')
            rep = menu('That File Alread Exists!','Replace the file','Save as a new file');
            if rep == 2
                while exist(filename,'file')
                    copy_num = copy_num + 1;
                    copy_txt = sprintf('(%d)',copy_num);
                    filename = [subFolder,'//','Subject_',Subject_str,'_T',trial_num_str,'_',copy_txt,num2str(now),'.mat'];
                    rawfilename = [subFolder,'//','Subject_',Subject_str,'_T',trial_num_str,'_',copy_txt,num2str(now),'_RAW.mat'];
                end
            end
        end
    end
    
    test_date = date;
    
    save(filename,'xpos','ypos','h','hdot','fuel','VN','VE','pitch','roll',...
        'crash','abort','distance2MeanSP', 'distance2Best_LP','LP_selection','xLP','yLP',...
        'reaction_times','rms_pitch','rms_roll','rms_total','Avail_Mode',...
        'trial_order','Subject','test_date','optimal','best_LP',...
        'differenceFromOptimal','num_selections','mean_rt','FPM','misses','alarm_time_1','Auditory_TIME_BUTTON')
%     % Save Audio File
%     audioFileName = strcat(filename,'audio.flac');
%     movefile ('AudioOutput.flac',audioFileName);

    clearvars -except rawfilename
    
    load('rawtmp.mat')

    save(rawfilename)
    
    delete('rawtmp.mat')
    
    clearvars -except Subject trial_num scenario_num use_ald hazards abort_press... 
    fail_detection fuel h hdot pitch roll Time VE VN xpos ypos valid_abort...
    train trials_done xLP yLP Avail_Mode thetaG phiG fuel_order trial_order...
    subFolder Subject_str safety LPs SPs RealTime LP_button IPC_button alarm_time_1 Auditory_TIME_BUTTON
else
    clear
    
    load('rawtmp.mat')
    
    delete('rawtmp.mat')
    
    clearvars -except Subject trial_num scenario_num use_ald hazards abort_press... 
    fail_detection fuel h hdot pitch roll Time VE VN xpos ypos valid_abort...
    train trials_done xLP yLP Avail_Mode thetaG phiG fuel_order trial_order...
    subFolder Subject_str safety LPs SPs RealTime LP_button IPC_button alarm_time_1 Auditory_TIME_BUTTON
end

safety = true;

