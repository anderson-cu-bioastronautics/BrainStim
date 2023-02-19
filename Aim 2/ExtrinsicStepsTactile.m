function [ON,Mag,dur_n]= ExtrinsicStepsTactile(u,dur,first_buzz_start_time,sec_time)
global stepnum second
ON=0;
Mag=0;
dur_n=dur; 
first_or_second=first_buzz_start_time;
if second==1
    first_or_second=sec_time;
end 

magnitudes=[.1:.025:.55];
times_on=first_or_second+linspace(0,length(magnitudes)*(dur),length(magnitudes));
times_on=times_on(1:length(magnitudes));
% if u  >= first_buzz_start_time-.1  && stepnum<=length(magnitudes)&&stepnum<=9
%     ON=1;
%     Mag=magnitudes(stepnum);
%     if stepnum==1
%         dur_n=dur+.1;
%     end
%     
% end
% end 
if (abs(u  - times_on(stepnum))<.1) && stepnum<length(magnitudes) && u<=times_on(end)
    ON=1;
    Mag=magnitudes(stepnum);
    if stepnum==1
        dur_n=dur+.1;
    end
elseif stepnum>=length(magnitudes) && second==0
    second=1;
    stepnum=1;
end 



end 