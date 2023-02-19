function []=TactileCall(stim_level,freq,duration)
string = sprintf('tactileHardware.py %0.7f %.7 %0.7f', abs(stim_level) , freq, duration );
command = string;
status = system(command);
end 