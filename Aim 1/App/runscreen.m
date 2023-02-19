function [outputint] = runscreen(stimlevel,screens,screenNumber,white,grey,black,inc,window,windowRect,app)
%negative is interval 1
%positive is interval 2

    display(stimlevel)
    display('interval 1 is negative, interval 2 is positive')
    display('Above is the stim level. This is an example scipt, this will not show in the actual script')
    

    
    
    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);

    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);

    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);

    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    % 0 is vertical, 90 is horizontal
    str={'0' '90'};
    if stimlevel < 0
        orientation = str2double(str(1));
    else
        orientation = str2double(str(2));
    end
    
% Play warning tone for user
% tedious and stupid, but needed if no noise is being played
if isempty(app.cha)==1
    app.cha=chami;
    app.cha.open();
    app.cha.abort_exams();
end

% empty and generate tone exam every single time
tg_exam=[];
tg_exam = tone_generation_exam();
ChannelOut = {[]  'HPR0'  []  []}; % Plays in right ear
Duration = 250;
tg_exam.ToneGenerationLevel.OutputChannel = ChannelOut;
tg_exam.ToneGenerationLevel.ToneDuration = Duration;
tg_exam.ToneGenerationLevel.Level = 80;
app.cha.queue_exam(tg_exam);
%     res = 22050;
%     len = 0.5*res;
%     hz = 440;
%     
%     sound(sin(hz*(2*pi*(0:len)/res)),res);
     pause(1)
    
    
    repeat = 1;
    while repeat ~= 3
        
    % Grating size in pixels
    gratingSizePix = 1000;

    % Grating frequency in cycles / pixel
    freqCyclesPerPix = .01;

    % Drift speed cycles per second
    cyclesPerSecond = 0;

    % Define Half-Size of the grating image.
    texsize = gratingSizePix / 2;

    % First we compute pixels per cycle rounded to the nearest pixel
    pixPerCycle = ceil(1 / freqCyclesPerPix);

    % Frequency in Radians
    freqRad = freqCyclesPerPix * 2 * pi;

    % This is the visible size of the grating
    visibleSize = 2 * texsize + 1;

    % Define our grating. Note it is only 1 pixel high. PTB will make it a full
    % grating upon drawing
    x = meshgrid(-texsize:texsize + pixPerCycle, 1);
    grating = grey * cos(freqRad*x) + grey;

    % Make a two layer mask filled with the background colour
    mask = ones(1, numel(x), 2) * grey;

    % Try to make horizontal
    mask1(:,:,1)=mask(:,:,2);
    mask1(:,:,2)=mask(:,:,1);
    mask=[];
    mask=mask1;

    % Contrast for our contrast modulation mask: 0 = mask has no effect, 1 = mask
    % will at its strongest part be completely opaque frameCounter.e. 0 and 100% contrast
    % respectively
    contrast = abs(stimlevel);

    % Place the grating in the 'alpha' channel of the mask
    mask(:, :, 2)= grating .* contrast;

    % Make our grating mask texture
    gratingMaskTex = Screen('MakeTexture', window, mask);

    % Make a black and white noise mask half the size of our grating. This will
    % be scaled upon drawing to make a "chunky" noise texture which our grating
    % will mask. Note the round function in here. For this demo we are simply
    % rounding the size to the nearest pixel, leaving PTB to do some scaling.
    noise = rand(round(visibleSize / 2)) .* white;

    % Make our noise texture
    noiseTex = Screen('MakeTexture', window, noise);

    % Make a destination rectangle for our textures and center this on the
    % screen
    dstRect = [0 0 visibleSize visibleSize];
    dstRect = CenterRect(dstRect, windowRect);

    % We set PTB to wait one frame before re-drawing
    waitframes = 1;

    % Calculate the wait duration
    waitDuration = waitframes * ifi;

    % Recompute pixPerCycle, this time without the ceil() operation from above.
    % Otherwise we will get wrong drift speed due to rounding errors
    pixPerCycle = 1 / freqCyclesPerPix;

    % Translate requested speed of the grating (in cycles per second) into
    % a shift value in "pixels per frame"
    shiftPerFrame = cyclesPerSecond * pixPerCycle * waitDuration;

    % Sync us to the vertical retrace
    vbl = Screen('Flip', window);

    % Set the frame counter to zero, we need this to 'drift' our grating
    frameCounter = 0;

    % Loop until a key is pressed
        % Calculate the xoffset for our window through which to sample our
        % grating
        xoffset = mod(frameCounter * shiftPerFrame, pixPerCycle);

        % Now increment the frame counter for the next loop
        frameCounter = frameCounter + 1;

        % Define our source rectangle for grating sampling
        srcRect = [xoffset 0 xoffset + visibleSize visibleSize];

        % Draw noise texture to the screen
        Screen('DrawTexture', window, noiseTex, [], dstRect, []);

        % Draw grating mask
        Screen('DrawTexture', window, gratingMaskTex, srcRect, dstRect, orientation);

        % Flip to the screen on the next vertical retrace
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

        pause(1);
        
        Screen('Flip', window);
        
        repeat = repeat + 1;
        if orientation == 0
            orientation = 90;
        else
            orientation = 0;
        end
        
        if repeat == 2
            pause(1);
        else
            break
        end
        
    end
        
      




                %this section just asks the current Matlab user for input

        %please replace this with a script that can adminster the stim in the
        %appropriate interval (if 'stimlevel' variable in negative, put stim in
        %interval 1, otherwise put it in interval 2
        %return subjects response (-1 for interval 1, +1 for interval 2)
        proper_response = 0;
        
        app.curStim=stimlevel;
        
        app.Lamp.Color=[0 1 0]; 
        
        if stimlevel < 0
            app.Presentation1Label.Visible = 1;
            app.Presentation2Label.Visible = 0;
            app.response=1;
        else
            app.Presentation1Label.Visible = 0;
            app.Presentation2Label.Visible = 1;
            app.response=2;
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
    proper_response = 1;
else
    % interval 2 response key
    outputint = 1;
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
            %app.confidence=50;
            saveData(app,app.WriteFile,app.n,app.curStim,app.response,app.correct(app.n));
            
%             while proper_response == 0
%                 %l_or_r_text = input('Response...Interval 1(1) or Interval 2(2)?:', 's');
%                 disp('Up is negative, right is positive')
%                 w = waitforbuttonpress;
%                 value = double(get(gcf,'CurrentCharacter'));
%                 if value == 30 %ID for the up button
%                     % interval 1 response key
%                     outputint = -1;
%                     proper_response = 1;
%                     loop = 1;
%                 elseif value == 29
%                     % interval 2 response key
%                     outputint = 1;
%                     proper_response = 1;
%                     loop = 1;
%                 else
%                     disp('WRONG KEY. PLEASE TRY AGAIN');
%                 end
%             end

      



end