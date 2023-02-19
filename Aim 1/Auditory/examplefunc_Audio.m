function [outputint] = examplefunc(stimlevel,tg_exam)
%negative is interval 1
%positive is interval 2

    display(stimlevel)
    OLevel=log10(abs(stimlevel));
    tg_exam.ToneGenerationLevel.Level = OLevel;
    if stimlevel < 0
        stim_interval = 1;
    else
        stim_interval = 2;
    end
    
    switch stim_interval
        case 1 % Play tone in first interval
            % Play tone
       app.Presentation1Label.BackgroundColor = [0 0 1];
       pause(.25/2);
       cha.queue_exam(tg_exam);
       pause(.25/2);
            % Pause for intermittant period
       pause(.25)
            % Disp second interval, without noise though
       disp(num2str(2));
       pause(.5)
       
        case 2 % Play tone in second interval
            % Disp second interval, without noise though
       disp(num2str(1));
       pause(.5);
            % Pause for intermittant period
       pause(.25)     
            % Play tone
       disp(num2str(2));
       pause(.25/2);
       cha.queue_exam(tg_exam);
       pause(.25/2);
    end
    %display('interval 1 is negative, interval 2 is positive')
    %display('Above is the stim level. This is an example scipt, this will not show in the actual script')
    
    
    %this section just asks the current Matlab user for input
    
    %please replace this with a script that can adminster the stim in the
    %appropriate interval (if 'stimlevel' variable in negative, put stim in
    %interval 1, otherwise put it in interval 2
    %return subjects response (-1 for interval 1, +1 for interval 2)
    proper_response = 0;
        while proper_response == 0
            l_or_r_text = input('Response...Interval 1(1) or Interval 2(2)?:', 's');
            if strcmp(l_or_r_text, '1')
                % interval 1 response key
                outputint = -1;
                proper_response = 1;
            elseif strcmp(l_or_r_text, '2')
                % interval 2 response key
                outputint = 1;
                proper_response = 1;
            else
                display('WRONG KEY. PLEASE TRY AGAIN');
            end
        end

end