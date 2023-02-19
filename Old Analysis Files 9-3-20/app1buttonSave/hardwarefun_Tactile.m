function [response] = hardwarefun_Tactile(stim_level,app)

set(app.DialogApp.RespondLabel,'visible','off')

app.curStim=stim_level;

%Visual warning
app.DialogApp.Presentation1Label.BackgroundColor = [0 0 1];
app.DialogApp.Presentation2Label.BackgroundColor = [0 0 1];
pause(.25);
app.DialogApp.Presentation1Label.BackgroundColor = [0.8 0.8 0.8];
app.DialogApp.Presentation2Label.BackgroundColor = [0.8 0.8 0.8];
pause(1);

if stim_level > 0
    
    app.DialogApp.Presentation1Label.BackgroundColor = [0 0 1];
    pause(0.1);
    % tactile stim in interval 1
    string = sprintf('/Users/BrainStim/Desktop/corbus-python/corbus/tactileHardware.py %0.7f %d', stim_level, 5);
    command = string;
    status = system(command);
    
    pause(0.2);
    app.DialogApp.Presentation1Label.BackgroundColor = [0.8 0.8 0.8];
    
    pause(1)
    
    app.DialogApp.Presentation2Label.BackgroundColor = [0 0 1];
    pause(0.1);
    string = sprintf('/Users/BrainStim/Desktop/corbus-python/corbus/tactileHardware.py %0.7f %d', 0.004, 5);
    command = string;
    status = system(command);
    
    pause(0.2);
 
    app.DialogApp.Presentation2Label.BackgroundColor = [0.8 0.8 0.8];
    pause(0.05)
else
    app.DialogApp.Presentation1Label.BackgroundColor = [0 0 1];
    pause(0.1);
    % tactile stim in interval 1
    string = sprintf('/Users/BrainStim/Desktop/corbus-python/corbus/tactileHardware.py %0.7f %d', 0.004, 5);
    command = string;
    status = system(command);
    
    pause(0.2);
    app.DialogApp.Presentation1Label.BackgroundColor = [0.8 0.8 0.8];
    
    pause(1)
    
    app.DialogApp.Presentation2Label.BackgroundColor = [0 0 1];
    pause(0.1);
    string = sprintf('/Users/BrainStim/Desktop/corbus-python/corbus/tactileHardware.py %0.7f %d', abs(stim_level), 5);
    command = string;
    status = system(command);
    
    pause(0.2);
 
    app.DialogApp.Presentation2Label.BackgroundColor = [0.8 0.8 0.8];
    pause(0.05);
end

set(app.DialogApp.RespondLabel,'visible','on') % Do something subject

%show correct answer on test admin's screen
if stim_level>0
    app.Presentation1Label.Visible=1;
    app.Presentation2Label.Visible=0;
else
    app.Presentation1Label.Visible=0;
    app.Presentation2Label.Visible=1;
end

%make lamp green to signify to admin that code is ready for subject's answer to be input
app.Lamp.Color=[0 1 0];

uiwait(); % Don't do anything until told to resume

app.Presentation1Label.Visible = 0;
app.Presentation2Label.Visible = 0;

if  app.Presentation1Button.Value == 1
    app.response = 1;
    response=1;
elseif app.Presentation2Button.Value == 1
    app.response = 2;
    response=-1;
else
    app.response=0;
    response=0;
end

if response*stim_level > 0
    % CORRECT!
    correct=1;
    app.correct(app.n)=correct;
    
    %concatenate previous stim level to prevStim vector
    app.prevStim(app.n)=app.curStim;
    
    %save new info
    if app.n==1
        app.time=0;
        app.start=tic;
    else
        app.time=toc(app.start);
        app.start=tic;
    end
    saveData(app,app.WriteFile,app.n,abs(app.curStim),app.response,app.correct,app.time)
    
    %update plot
    plotpts=abs(app.prevStim);
    stairs(app.UIAxes,[1:app.n],plotpts,'-o');
    
elseif response*stim_level < 0
    correct=0;
    app.correct(app.n)=correct;
    
    %concatenate previous stim level to prevStim vector
    app.prevStim(app.n)=app.curStim;
    
    %save new info
    if app.n==1
        app.time=0;
        app.start=tic;
    else
        app.time=toc(app.start);
        app.start=tic;
    end
    saveData(app,app.WriteFile,app.n,abs(app.curStim),app.response,app.correct,app.time)
    
    %update plot
    plotpts=abs(app.prevStim);
    stairs(app.UIAxes,[1:app.n],plotpts,'-o');

end

%make lamp yellow again for next trial to start
app.Lamp.Color=[1 1 0];

uiwait();
%change lamp to red to signify the stims are being presented
app.Lamp.Color=[1 0 0];
%reset the answer so the lapse button is selected
app.LapseButton.Value=1;


end
