function [outputint] = noappscreen(stimlevel,screens,screenNumber,white,grey,black,inc,window,windowRect)
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
        c = clock;
        initTime = clock;
        newTime = initTime;
        while newTime(6)-initTime(6) < 1
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
          newTime = clock;
        end
        
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
            
            while proper_response == 0
                %l_or_r_text = input('Response...Interval 1(1) or Interval 2(2)?:', 's');
                disp('Up is negative, right is positive')
                w = waitforbuttonpress;
                value = double(get(gcf,'CurrentCharacter'));
                if value == 30 %ID for the up button 
                    % interval 1 response key
                    % Vertical = negative response
                    % Vertical was first
                    outputint = -1;
                    proper_response = 1;
                    loop = 1;
                elseif value == 29
                    % interval 2 response key
                    % Horizontal = positive response
                    % Horizontal was first
                    outputint = 1;
                    proper_response = 1;
                    loop = 1;
                else
                    disp('WRONG KEY. PLEASE TRY AGAIN');
                end
            end

      



end