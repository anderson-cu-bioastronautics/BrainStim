function out=usePython(time,tactor,Amp,dur)
% This function is used to call CallTactile which runs a python script 
% 'tactileHardware.py' and establishes serial communication with the piezo
% tactile device. 
% This function also defines extrinsic global variables to don't repeat the
% same buzz twice and do not call the tactor device twice in a given
% instance.    WHEN MODIFIYING CODE open the simulink diagnostics and DO
% NOT CONNECT THE TACTOR, and instead look for a message "Could not find any CorBus USB interfaces! "
% which is equivalent to a buzz when the tactor is disconnected. 


% Daniel Gutierrez-Mendoza

out=0;
global extr myExtrinsicTime stepnum
if time>myExtrinsicTime
    extr=0;
end 
if tactor==1 && extr==0 
    % If tactor is on then the thing will command it to buzz. 
    CallTactile(Amp,5,dur);
%     extr=1;
%     myExtrinsicTime=time;
    stepnum=stepnum+1;
  
end 


end 