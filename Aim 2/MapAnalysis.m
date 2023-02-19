% First load your map from the workspace, so right click the map's .mat
% file and hit load

% Calculate the mean science point
lander_size = 50; %[ft]
numPoints = 2000;
feetPerPoint = 24000/numPoints;
SP_mean = [24000 - mean(SPs(2,:)) ; 24000 - mean(SPs(1,:))]; % calculate SP centroid

% Calculate the distance for each landing point
for n = 1:6
    LP = [LPs(1,n);24000 - LPs(2,n)]; %[feet]
    dist_SP_NH(n) = norm(SP_mean - LP); %[feet] distance of LP from centroid
%    dist_LP(n) = norm(LP - [Xffeet;24000-Yffeet]); %[feet] not needed
    
    % Eliminate Hazardous Points
    LP = LP/feetPerPoint; %[points]
    X_LP = (LP(1) - round(lander_size/(2*feetPerPoint)):(LP(1) + round(lander_size/(2*feetPerPoint)) - 1)); %[points]
    Y_LP = (2000 - LP(2) - round(lander_size/(2*feetPerPoint))):(2000 - LP(2) + round(lander_size/(2*feetPerPoint)) - 1); %[points]
    
    LP_hazards = 100*sum(sum(hazards(Y_LP,round(X_LP))))/(lander_size/feetPerPoint)^2; %[percent]
    
    choice_SP(1,n) = n; % LP name... I think
    choice_SP(2,n) = dist_SP_NH(n); % distance from centroid
    if LP_hazards > 25
    choice_SP(3,n) = 1; % yes it is hazardous
    else
    choice_SP(3,n) = 0; % no it is not hazardous... lucky you
    end
end

% Rank the landing zones, post hazard. THIS DOES NOT RANK PRE-HAZARD! You
% sadly have to do this on your own
hazLP = choice_SP(:,find(choice_SP(3,:)==1));
goodLP = choice_SP(:,find(choice_SP(3,:)==0));
[~,sort_haz] = sort(hazLP(2,:));
[~,sort_good] = sort(goodLP(2,:));

rank_LP = [goodLP(:,sort_good) hazLP(:,sort_haz)]