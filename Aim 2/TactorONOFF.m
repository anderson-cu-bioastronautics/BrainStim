function [ON] = TactorONOFF(time)
ON=1;
global TactorTimes
if isempty(find(TactorTimes==time))==1
    ON=0;
end 
end

