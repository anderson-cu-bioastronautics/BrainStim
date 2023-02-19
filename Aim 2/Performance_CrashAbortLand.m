%% Crash, Abort, Land Performance Metric
clc; 

Subject = 1000;
if Subject >9
    Subject_str = num2str(Subject);
else
    Subject_str = ['0',num2str(Subject)];
end
subFolder = ['Subject_Data//Subject_',Subject_str];

trial = 4;

if trial >9
        trial_str = num2str(trial);
    else
        trial_str = ['0',num2str(trial)];
end
file_str = ['Subject_' Subject_str '_T' trial_str '_RAW.mat'];
filename = fullfile([pwd '//' subFolder '//' file_str]);
load(filename)


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

fprintf('Crash, Abort, Land Score: %.3f \n\n',CAL_Score)

