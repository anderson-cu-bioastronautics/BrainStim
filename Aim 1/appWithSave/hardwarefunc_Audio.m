% script to randomize stimulus presentation

function [outputint] = hardwarefunc_Audio(stimlevel,app)

% tedious and stupid, but needed if no noise is being played
if isempty(app.cha)==1
    app.cha=chami;
    app.cha.open();
    app.cha.abort_exams();
end

%negative is interval 1
%positive is interval 2
set(app.DialogApp.RespondLabel,'visible','off')
%    display(stimlevel)

% empty and generate tone exam every single time
tg_exam=[];
tg_exam = tone_generation_exam();
ChannelOut = {[]  'HPR0'  []  []}; % Plays in right ear
Duration = 250;
tg_exam.ToneGenerationLevel.OutputChannel = ChannelOut;
tg_exam.ToneGenerationLevel.ToneDuration = Duration;
OLevel=log10(abs(stimlevel));
tg_exam.ToneGenerationLevel.Level = OLevel;
app.curStim=OLevel;
app.StimDisplayLabel.Text=[num2str(OLevel) ' dB'];

if stimlevel < 0
    stim_interval = 1;
else
    stim_interval = 2;
end

switch stim_interval
    case 1 % Play tone in first interval
        % Play tone
        app.DialogApp.Presentation1Label.BackgroundColor = [0 0 1];
        pause(.25/2);
        app.cha.queue_exam(tg_exam);
        pause(.25/2);
        app.DialogApp.Presentation1Label.BackgroundColor = [0.8 0.8 0.8];
        % Pause for intermittant period
        pause(.25)
        % Disp second interval, without noise though
        app.DialogApp.Presentation2Label.BackgroundColor = [0 0 1];
        pause(.5)
        app.DialogApp.Presentation2Label.BackgroundColor = [0.8 0.8 0.8];
    case 2 % Play tone in second interval
        % Disp first interval, without noise though
        app.DialogApp.Presentation1Label.BackgroundColor = [0 0 1];
        pause(.5)
        app.DialogApp.Presentation1Label.BackgroundColor = [0.8 0.8 0.8];
        % Pause for intermittant period
        pause(.25)
        % Play tone
        app.DialogApp.Presentation2Label.BackgroundColor = [0 0 1];
        pause(.25/2);
        app.cha.queue_exam(tg_exam);
        pause(.25/2);
        app.DialogApp.Presentation2Label.BackgroundColor = [0.8 0.8 0.8];
end

set(app.DialogApp.RespondLabel,'visible','on') % Do something subject


%this section just asks the current Matlab user for input

%please replace this with a script that can adminster the stim in the
%appropriate interval (if 'stimlevel' variable in negative, put stim in
%interval 1, otherwise put it in interval 2
%return subjects response (-1 for interval 1, +1 for interval 2)
proper_response = 0;
%make lamp green to signify to admin that code is ready for subject's answer to be input
app.Lamp.Color=[0 1 0];
%show correct answer on test admin's screen
            if stimlevel<0
                    app.Presentation1Label.Visible=1;
                    app.Presentation2Label.Visible=0;
            else
                    app.Presentation1Label.Visible=0;
                    app.Presentation2Label.Visible=1;
            end
uiwait(); % Don't do anything until told to resume
% app.state=0;
% while app.state == 0
%     waiting=1;
% end
% l_or_r_text = input('Response...Interval 1(1) or Interval 2(2)?:', 's');
if  app.Presentation1Button.Value == 1
    % interval 1 response key
    outputint = -1;
    app.response=1;
    proper_response = 1;
else
    % interval 2 response key
    outputint = 1;
    app.response=2;
    proper_response = 1;
    %             else
    %                 display('WRONG KEY. PLEASE TRY AGAIN');
    
end
        if outputint*stimlevel > 0
            % CORRECT!
            correct=1;
        else
            correct=0;
        end
            app.correct(app.n)=correct;
            app.confidence=50;
end

