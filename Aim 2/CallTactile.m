function out=CallTactile(Amp,Freq,Dur)
out=0;
%global stop;
if Amp<1 %&& stop ==0
    %string = sprintf('tactileHardware.py %0.7f %d %0.7f', abs(Amp) , Freq, Dur );
    %string = sprintf('USMESSING.py %0.7f %d %0.7f', abs(Amp) , Freq, Dur );
    system(string);
    out=1;
   % stop=1;
end 
end 