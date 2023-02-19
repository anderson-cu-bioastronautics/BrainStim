function [Xdot] = trajectory_sim(t, X) % CHANGE

% Modified: 4/18/17 - removed global variables

global R0val Rchecker Xdot other_outputs
% 
vars = ones(15,1);  % CHANGE
tolv = ones(7,1);   % CHANGE
hdotTD = -3;    % ft/s
hLG = 500;      % ft
hdotLG = -16;   % ft/s
hTD = 150;      % ft
y0 = 12*2000/2; 
yLP0 = y0 + 1350; 
x0 = 12*2000/2; 
xLP0 = x0;

R0 = sqrt((y0-yLP0)^2 + (x0-xLP0)^2);%vars(1);
hTD= 150;
tauV= 8;
tauh= 25;
gBody= 1.622 * 3.28084; 
tauT= 1.5;
Isp= 311;
g0= 9.81 * 3.28084;
a = (hdotLG - hdotTD)/(hLG-hTD);
b = -a*hTD + hdotTD;
C = -b/a;
k = hLG-C;
q = (2*a)/log(k/(hTD-C));
K = 1;
% the limit for how small the range has to be for us to continue to descend if we are below 150 ft
Rlim = 50;
R_R0_tol = exp(1)^2 - 1.01;
max_thrust = 10000; %lbf, defined by Bilimoria
thrust_u = .6;  % 60%, defined by Bilimoria
thrust_l = .1;  % 10%, defined by Bilimoria
max_phi = 45;   % deg, defined by Bilimoria
max_theta = 45; % deg, defined by Bilimoria


% Decompose the state vector
x = X(1);       % ft, horizontal North position (North of the LP is positive)
y = X(2);       % ft, horizontal East position (East of the LP is positive) 
VN = X(3);      % ft/s, horizontal North velocity (headed North is positive)
VE = X(4);      % ft/s, horizontal East velocity (headed East is positive)
theta = X(5);   % degrees, pitch angle (positive is nose up)
phi = X(6);     % degrees, roll angle (positive is right)
h = X(7);       % ft, height (above the surface is positive), including terrain at this point
hdot = X(8);    % ft/s, descent rate (up is positive, down is negative)
m = X(9);       % slugs, mass of the vehicle as fuel is burned


%% Perform simulation as normal
% Other states assumed to be 0
psi = 0;    % assume there is no yaw

R = sqrt(y^2+x^2);      % ft, range to LP, always positive, assume it is only one horizontal dimension (y=0)
V = sqrt(VN^2+VE^2);    % ft/s, net horizontal velocity, assume it is only one horizontal dimension

% Our added logic to deal with cases where you are farther away from LP
% than when you started the landing
if (R/R0) > (exp(1)^2 - R_R0_tol)
    % Change R0 to your actual range because VN blows up  if R/R0 is
    % greater than e^2. This will create a "new" simulation (new initial
    % range that you are evaluating from). Essentially  restarting the
    % simulation where the previous state left of
    % Set R0 global to the current range value (redefine R0) (new R0 for
    % the remaining time steps)
    R0val = R;
    % Chnages R0 for the current time step
    R0 = R;
    Rchecker = 1;
    
end

% to determine if range is increasing or decreasing, we need to decide if
% the horizontal velocity is headed towards the LP or away from it
angleLP = 180/pi*atan2(-x,-y);
Rdot = -VN*cosd(angleLP) - VE*sind(angleLP);  % ft/s, time rate of change in range (could be negative if moving away from LP)

% Hdot Guidance
if h >= hTD
    hG = C + (((hTD-C)^log(sqrt(R/R0)))/k)^(1/(log(sqrt(R/R0))-1));     % ft, Equation 10
    if hG > 550
        hG = 550;
    end
    hdotG = (a*hG + b) + (hG - h)/tauh;                               % ft/s, Equation 11
else
    hG = h;     % ft, not a real thing, just a place holder for once we get below 150 ft since hG is no longer defined
    hdotG = -3;                                                       % ft, if below 150 ft, just descend at -3 ft/s
end

% Add in a limit that if we are too far away (>Rlim) but too low (<150ft),
% go into hover (this is our own personal addition, not from Bilimoria)
if h < hTD && R > Rlim
    hdotG = 0;  % go into a hover (i.e. hdotG=0)
end

% Limit hdotG to between 0 and -30 f/s
if hdotG > 0
    hdotG = 0;
elseif hdotG < -30
    hdotG = -30;
end

% Define a deadband that will be added in on the Tdetlacmd when assum
hdot_deadband = 0.1;                % ft/s
% The commanded vertical velocity will always be the guidance vertical velocity
hdotcmd = hdotG;        % assuming this to be the case, need to verify
% thrust command
T0cmd = m*gBody/(cosd(theta)*cosd(phi));       % ft-lbs (like N), Equation 1, primary thrust command to counteract weight of the vehicle
% change in thrust command
Tdeltacmd = m/(cosd(theta)*cosd(phi)) * (hdotcmd-hdot)/tauT;    % ft-lbs (like N), Equation 2, secondary thrust command to adjust descent rate
% Add a deadband in on the Tdetlacmd
if abs(hdotcmd-hdot) < hdot_deadband
    Tdeltacmd = 0;
end
T = T0cmd + Tdeltacmd;                          % ft-lbs (like N), total thrust command
% Constrain thrust between 10% and 60% of max operational thrust (10,000 lbf)
if T > max_thrust*thrust_u
    T = max_thrust*thrust_u;
end
if T < max_thrust*thrust_l
    T = max_thrust*thrust_l;
end
% hdot = hdot;    % ft/s, duh.
hddot = T/m * cosd(phi)*cosd(theta) - gBody;    % ft/s^2, vertical acceleration
% Change in vehicle mass
mfueldot = T/(Isp*g0);      % slugs/s, fuel consumption rate, assumes ideal rocket equation
mdot = -mfueldot;           % slug/s

% Velocity guidance (North and East)
VNG = y*q*(1-log(sqrt(R/R0)));     % ft/s, Equation 13a
VEG = x*q*(1-log(sqrt(R/R0)));     % ft/s, Equation 13b

% acceleration Guidance (North and East)
% Set x acceleration equal to zero if x = 0
if y == 0
    aNG = 0;
else
    aNG = (VN*VNG/y + q*(-1)*y*Rdot/(2*R)) + (VNG - VN)/tauV;  % ft/s^2, Equation 15a
end
% Set y acceleration equal to zero if y = 0
if x == 0
    aEG = 0;
else
    aEG = (VE*VEG/x + q*(-1)*x*Rdot/(2*R)) + (VEG - VE)/tauV;  % ft/s^2, Equation 15b
end

% roll Guidance
phiG = asind(-m/T * (aNG*sind(psi) - aEG*cosd(psi)) );   % deg, (positive is nose up), Equation 17a
% Limit phiG by +/- 45 deg
if phiG > max_phi
    phiG = max_phi;
elseif phiG < -max_phi
    phiG = -max_phi;
end

% pitch Guidance
thetaG = asind(-m/(T*cosd(phiG)) * (aNG*cosd(psi) + aEG*sind(psi)) );   % deg, (positive is nose up), Equation 17b
% Limit thetaG by +/- 45 deg
if thetaG > max_theta
    thetaG = max_theta;
elseif thetaG < -max_theta
    thetaG = -max_theta;
end

% Set the derivatives
xdot = VN;      % ft/s, positive for moving North
ydot = VE;      % ft/s, positive for moving East

aN = -T/m*(cosd(phi)*sind(theta)*cosd(psi) + sind(phi)*sind(psi));  % ft/s^2, Equation 16a, positive for accelerating towards LP
VNdot = aN;           % ft/s^2
aE = -T/m*(cosd(phi)*sind(theta)*sind(psi) - sind(phi)*cosd(psi));  % ft/s^2, Equation 16b, positive for accelerating towards east
VEdot = aE;           % ft/s^2

% Pilot makes control inputs to align theta with thetaG similar to a
% proportional control with gain Ktheta. We assume no time delay.
% K = 0.5;     % proportional nulling gain. Currently selected arbitarily, need to updated this with actual limits and/or better info on pilot behavior.
thetadot = K * (thetaG - theta);   % deg/s

% Limit thetadot based upon how much control authority exists
thetadot_lim = 5;       % deg/s. Limit on how quickly theta can change per unit second. Currently selected arbitarily. Prevents large errors in pitch - pitch guidance to cause massive pitch rates which would not be instanteously relizable.
if thetadot > thetadot_lim
    thetadot = thetadot_lim;
elseif thetadot < -thetadot_lim
    thetadot = -thetadot_lim;
end

% Pilot makes control inputs to align theta with phiG similar to a
% proportional control with gain Ktheta. We assume no time delay.
% K = 0.5;     % proportional nulling gain. Currently selected arbitarily, need to updated this with actual limits and/or better info on pilot behavior.
phidot = K * (phiG - phi);   % deg/s

% Limit thetadot based upon how much control authority exists
phidot_lim = 5;       % deg/s. Limit on how quickly theta can change per unit second. Currently selected arbitarily. Prevents large errors in pitch - pitch guidance to cause massive pitch rates which would not be instanteously relizable.
if phidot > phidot_lim
    phidot = phidot_lim;
elseif phidot < -phidot_lim
    phidot = -phidot_lim;
end

% Accumulate other outputs
other_outputs = [hG; hdotG; T0cmd; Tdeltacmd; T; VNG; aNG; thetaG; VEG; aEG; phiG; mdot; thetadot; phidot];

% Recompose the derivative of the state vector
Xdot = [xdot; ydot; VNdot; VEdot; thetadot; phidot; hdot; hddot; mdot];