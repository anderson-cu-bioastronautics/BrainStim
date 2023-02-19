stim_level = 0.4;
string = sprintf('tactileHardware.py %0.7f %d', abs(stim_level) , 400);
command = string;
status = system(command);