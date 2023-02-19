function [] = hardwarefun_Tactile(stim_level,app)

clear a
a = arduino();

sound = 1; % switch for sound                        

% Buzzer set up 
buzzer_pin = 'D10'; 
sound_frequency = 262;

% LED set up
led = 'D6';

% Forced choice input set up 
configurePin(a, 'D4', 'pullup');
configurePin(a, 'D2', 'pullup');
button_1 = 'D4';
button_2 = 'D2';

%initalize other varibles and flags
stop = 0;                           % Amplitude, start with no vibrations
dac = i2cdev(a,'0x60');              % Connect to mcp4725 using I2C protocol
write(dac,stop,'uint16');

correct_interval = randi(2);

switch correct_interval 
    
    case 1
        app.DialogApp.Presentation1Label.BackgroundColor = [0 0 1];
        % light and sound
        if sound == 1
            playTone(a, buzzer_pin, sound_frequency, 0.2);
        end
        writeDigitalPin(a, led, 1);
        pause(0.1);
        writeDigitalPin(a, led, 0);                            

        pause(0.15);                                              
        
        % tactile stim in interval 1
        write(dac,stim_level,'uint16')                          
        pause(1.2)
        write(dac,stop,'uint16')
        %fprintf('Stim Administered: %.2f \n', stim_level); 

        pause(0.03);
        
        app.DialogApp.Presentation1Label.BackgroundColor = [0.8 0.8 0.8];
        app.DialogApp.Presentation2Label.BackgroundColor = [0 0 1];
        
        % light and sound 
        if sound == 1
            playTone(a, buzzer_pin, sound_frequency, 0.2);
        end  
        writeDigitalPin(a, led, 1);
        pause(0.1);
        writeDigitalPin(a, led, 0);                            
        
        pause(0.15);                                            
        
        % no tactile stim in interval 2
        pause(1.2);                                             
        
        if sound == 1
            playTone(a, buzzer_pin, sound_frequency, 0.1);
            pause(0.1)
            playTone(a, buzzer_pin, sound_frequency, 0.1);
        end
        
        app.DialogApp.Presentation2Label.BackgroundColor = [0.8 0.8 0.8];
        
        
    case 2
        app.DialogApp.Presentation1Label.BackgroundColor = [0 0 1];
        % light and sound                                       
        if sound == 1
            playTone(a, buzzer_pin, sound_frequency, 0.2);
        end
        writeDigitalPin(a, led, 1);
        pause(0.1);
        writeDigitalPin(a, led, 0);                             
        
        pause(0.15);                                               
        
        % no tactile stim in interval 1
        pause(1.2);                                              
        
        app.DialogApp.Presentation1Label.BackgroundColor = [0.8 0.8 0.8];
        app.DialogApp.Presentation2Label.BackgroundColor = [0 0 1];
        % light and sound 
        if sound == 1
            playTone(a, buzzer_pin, sound_frequency, 0.2);
        end 
        writeDigitalPin(a, led, 1);
        pause(0.1);
        writeDigitalPin(a, led, 0);                          
        
        % tactile stim in interval 2
        write(dac,stim_level,'uint16')                         
        pause(1.2)
        write(dac,stop,'uint16');
        %fprintf('Stim Administered: %.2f \n', stim_level); 
        
        if sound == 1
            playTone(a, buzzer_pin, sound_frequency, 0.1);
            pause(0.1)
            playTone(a, buzzer_pin, sound_frequency, 0.1);
        end
       
        app.DialogApp.Presentation2Label.BackgroundColor = [0.8 0.8 0.8];
        set(app.DialogApp.RespondLabel,'visible','on') % Do something subject

%make lamp green to signify to admin that code is ready for subject's answer to be input
app.Lamp.Color=[0 1 0]; 
%show correct answer on test admin's screen
            switch correct_interval
                case 1
                    app.Presentation1Label.Visible=1;
                    app.Presentation2Label.Visible=0;
                case 2
                    app.Presentation1Label.Visible=0;
                    app.Presentation2Label.Visible=1;
            end
uiwait(); % Don't do anything until told to resume


if  app.Presentation1Button.Value == correct_interval
    correct = 1;
else
    correct = 0;
end
    app.correct(app.n)=correct;
    app.confidence=50;
    saveData(app,app.filename,app.n,app.curStim,app.correct(app.n),app.confidence);

end
