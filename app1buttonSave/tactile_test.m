stim_level = 0.04;
string = sprintf('/Users/BrainStim/Desktop/Tactor_Hardware/corbus-python-master/tactileHardware.py %0.7f %d', abs(stim_level) , 400);
command = string;
status = system(command);