function [value,isterminal,direction] = stopping_point(tOUT, X)
% Locate the time when height passes through zero in a decreasing direction and stop integration.

% detect when landed
value(1) = X(7); % detect when h = X(7) < 0
isterminal(1) = 1; % stop the integration
direction(1) = -1; % negative direction

value(2) = X(9) - 493; % detect when h = X(7) < 0
isterminal(2) = 1; % stop the integration
direction(2) = -1; % negative direction


