function mepFunction
clear all
%Channel list should be
%Channel 0- TMS
%Channel 1- FDI (target muscle for MEPs)
%Channel 2- Muscle 2
%Channel 3- Muscle 3
%Channel 4- Muscle 4
clear all;
delete(instrfindall);
ListenChar(0); %to ensure that the keyboard is responsive
backgnd_EMG_threshold=0.007; %Above this it cant move from state 1 to 2 or 2 to 3
% % if using matlab and psychtoolbox is installed, set to highest priority
% % (9 for mac; 2 for windows)
 Priority(2);
% % initialise virtual com (serial) port
 usb = serial('COM3'); % virtual com port address will look something like this on windows
% % set necessary parameters; any program you use (eg presentation) should
% allow you to specifiy these options
% set(usb,'BaudRate',115200);
% set(usb,'DataBits',8);
% set(usb,'Parity','none');
% set(usb,'StopBits',1);
% marker = 1; % the marker can be any value between 1 and 254
% % open serial port
% fopen(usb);
% first write is always slow, so perform a dummy write before experiment
% starts
% fwrite(usb,1); 
%%
nMEPs = 30;  % number of MEPs to collect

rmt= input('Desired Intensity');
% set experiment parameters
subNum= input('Subject number');
subCode = ['NF' num2str(subNum)];
direc = 'C:\Users\big lab\Documents\MATLAB';
blockNumber = input('Block Number?');
if blockNumber>4
    FirstOrSecondDay=2;
else
    FirstOrSecondDay=1;
end
mean_MEP=input('Baseline MEP size?');
min_MEP_value=0; max_MEP_value=mean_MEP*2; Range=max_MEP_value-min_MEP_value;
% Seed random number generator
rng('shuffle');
%load('C:\Users\big lab\Documents\MATLAB\MEP\dummy_rms_data');
%In control condition, generate vector of MEPs to read in. Frequency of
%reward is average from previous experiment
percentage_reward=0.5; % 50% of trials will have a reward
%number_rewards=nMEPs*percentage_reward;
number_rewards=20;
number_nonrewards=nMEPs-number_rewards;
bigger_than_mean_fake_MEPS=((mean_MEP*2)-mean_MEP).*rand(number_rewards,1)+mean_MEP; %swap the rewards and nonrewards variables in UP and DOWN scripts
smaller_than_mean_fake_MEPS=(mean_MEP-(mean_MEP*2)).*rand(number_nonrewards,1)+mean_MEP;

control_MEP_vector=zeros(nMEPs,1);

%% I Think I need to find out exactly what the reward percentage will be and build a vector of 1 or 0s for it. then can just change this
%%assuming reward frequency is 50%
reward_vector=[0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1];
reward_positions=find(reward_vector==1)
no_reward_positions=find(reward_vector==0)
%load('C:\Users\big lab\Documents\MATLAB\MEP\dummy_rms_data');
b=1;
for pos=1:length(reward_positions)
    control_MEP_vector(reward_positions(pos))=bigger_than_mean_fake_MEPS(b)  %specific to UP condition
    b=b+1
end
b=1;
for pos=1:length(no_reward_positions)
    control_MEP_vector(no_reward_positions(pos))=smaller_than_mean_fake_MEPS(b) % specific to UP condition
    b=b+1
end    

nChans = 5;  % number of EMG channels
emgDuration = 1;  % collect EMG data for this duration
triggerTime = 0.5;  % time after trial start to send TMS trigger
feedbackDuration = 4;  % visual feedback duration
minGap = 3; maxGap = 6;   % set min and max gap times
gapDurationArray = (rand(nMEPs,1)*maxGap) + minGap;  % array of gap duration between trials
minfixGap = 1; maxfixGap = 3;   % set min and max gap times
fixgapDurationArray = (rand(nMEPs,1)*maxfixGap) + minfixGap;  % array of gap duration between trials


% file handling
if FirstOrSecondDay ==1
    SesNum='Ses1';
    subDir = [direc filesep 'data' filesep subCode filesep 'Ses1'];
    if ~exist(subDir,'dir')
        mkdir(subDir);
    end
    
else
    SesNum='Ses2';
    subDir = [direc filesep 'data' filesep subCode filesep 'Ses2'];
    if ~exist(subDir,'dir')
        mkdir(subDir);
    end
end

dataFilename = [direc filesep 'data' filesep subCode filesep SesNum filesep 'UP_b'  num2str(blockNumber) '-m' num2str(nMEPs) '-c' num2str(nChans) '.dat'];

%This is protection for overwriting files. Keep commented for programming
%purposes but uncomment during real experiment
    if fopen(dataFilename,'rt') ~= -1
        fclose('all');
        error('computer says no: result data file already exists!');
    else
    dataFilePointer = fopen(dataFilename,'wt');
    end

% define global variables
global rD
global tS

%%% configure TMS serial port communication here %%%

% setup data acquisition device
s = daq.createSession('ni');
s.addAnalogInputChannel('Dev2', 0:(nChans-1), 'Voltage');
set(s.Channels, 'InputType', 'SingleEnded');
set(s.Channels, 'Range', [-5,5]);
s.Rate = 2000;
s.NotifyWhenDataAvailableExceeds = 60;
s.IsContinuous = true;
lh = addlistener(s, 'DataAvailable', @rData);
s.startBackground();

% Call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Give Matlab high priority
Priority(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% The avaliable keys to press
escapeKey = KbName('ESCAPE');
spaceKey = KbName('SPACE');

% Get the centre coordinate of the window
[centreX, centreY] = RectCenter(windowRect);

% Enable alpha blending for anti-aliasing
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Hide cursor
HideCursor;

% Sync us and get a time stamp
vbl = Screen('Flip', window);
waitframes = 1;

% Maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

% Dummy calls to make sure functions are ready to go without delay
KbCheck;
[keyIsDown, secs, keyCode] = KbCheck;
WaitSecs(0.1);
GetSecs;

% Set default screen font and size for written messages
Screen('TextSize', window, 50);
Screen('TextFont', window, 'Arial');
dotSizePix=80;

%% SET UP TMS AND SERIAL PORT COMMUNICATION
%Set up serial port connection
serialPortObj=serial('COM1', 'BaudRate',9600,'DataBits',8,'Stopbits',1,'Parity','none','FlowControl','none','inputbuffersize',1024,'outputbuffersize',1024,'Terminator','?');
% Callback function to execute every 500 ms to ensure that the stimulator
% is in the remote control mode and will stay armed. Otherwise,
% stimulator will disarm itself automatically in about 1 sec.
serialPortObj.TimerPeriod = 0.5; % period of executing the callback function in sec
fopen(serialPortObj);
serialPortObj.TimerFcn = {'Rapid2_MaintainCommunication'};
Rapid2_Delay(1000, serialPortObj);
pause on;
%arm stimulator
%success = Rapid2_ArmStimulator(serialPortObj);
%set power level
powerLevel=rmt;

if powerLevel>100;
    powerLevel=100;
end

success = Rapid2_SetPowerLevel(serialPortObj, powerLevel, 1);
if ~success
    display 'Error: Cannot set the power level';
    return
else
    % Display power level;
    display(powerLevel);
    % Introduce delay to allow the  stimulator to adjust to the new power level
    %Rapid2_Delay(4000, serialPortObj);
end
%%Prepare the sounds
% load sounds
soundsPath= 'C:\Users\big lab\Documents\MATLAB\MEP\sounds';
winWavPath = [soundsPath filesep 'better.wav'];
winSoundData = audioread(winWavPath);
winSoundData=winSoundData';

loseWavPath = [soundsPath filesep 'worse.wav'];
loseSoundData = audioread(loseWavPath);
loseSoundData = loseSoundData';
% Perform basic initialization of the sound driver:
InitializePsychSound(1);
pahandle = PsychPortAudio('Open', 3, [], 2,[], 2, 0, 0);

%% Create a line that will be drawn at the mean MEP size
%Calculate the range of MEP sizes
barColour=[0 1 0];
%Format is left, top, right, bottom
lineYpos=screenYpixels-round(screenYpixels*(mean_MEP/Range));
line=[0 lineYpos screenXpixels lineYpos+10];
lineColour=[1 1 1];
MEP_amp=0;

%Add textures
tick_path= 'C:\Users\big lab\Documents\MATLAB\MEP\tick_transparent_smaller.png';
[tick map tickA] = imread(tick_path);
tick(:,:,4) = tickA;
[tickX tickY tickD] = size(tick);
tickTexture = Screen('MakeTexture', window, tick);
cross_path= 'C:\Users\big lab\Documents\MATLAB\MEP\cross2_small.png';
[cross map crossA] = imread(cross_path);
cross(:,:,4) = crossA;
[crossX crossY crossD] = size(cross);
crossTexture = Screen('MakeTexture', window, cross);
targetX=200;
targetY=200;
dollar_path= 'C:\Users\big lab\Documents\MATLAB\MEP\Doublebarred_dollar_sign.png';
[dollar map dollarA] = imread(dollar_path);
dollar(:,:,4) = dollarA;
[dollarX dollarY dollarD] = size(dollar);
dollarTexture = Screen('MakeTexture', window, dollar);
dollar_targetX=centreX+300;
dollar_targetY=900;
fixation_path= 'C:\Users\big lab\Documents\MATLAB\MEP\white_plus_smaller.png';
[fixation map fixationA] = imread(fixation_path);
fixation(:,:,4) = fixationA;
[fixationX fixationY fixationD] = size(fixation);
fixationTexture = Screen('MakeTexture', window, fixation);
fixation_targetX=centreX;
fixation_targetY=centreY;
greenDot=[0 1 0];
%% START EXPERIMENT

DrawFormattedText(window, 'press the spacebar to begin experiment', 'center', 'center', white);
Screen('Flip', window);
while (keyCode(spaceKey) == 0) [keyIsDown, secs, keyCode] = KbCheck; end
keyCode(spaceKey) = 0;
success = Rapid2_ArmStimulator(serialPortObj);
WaitSecs(0.5);
nSuccess=0;
mepCount = 0;

rawRMSData = [];
dummyData=zeros(1000,6); %1000 samples is 500ms
dummyData(:,:)=0.005;
dummyData(:,1)=rand(1,length(dummyData)); %I added this becasue when I later use the Unique function for the first time on the raw data, if the timestamps are all identical, it cancels out all data as they are not unique. This is only needed for the dummy samples that initialise the matrix that later gets filled with rms data
rawRMSData=[rawRMSData; dummyData];
    FDIrmsgood=1;
    ch3rmsgood=1;
    ch4rmsgood=1;
    ch5rmsgood=1;
    
    %Here I added code to collect a small chunk of data to estimate the
    %bias that needs to be removed from each EMG channel. Very important to
    %do this for accurate RMS measurements
    bias_estimation_duration=2;
    bias_starttime=0;
    rawchunkData =[];
    while bias_starttime<bias_estimation_duration
        bias_starttime = GetSecs;
        timeStamps = tS;
        rawData = rD;
         rawchunkData = [rawchunkData; timeStamps rawData];
    end
    uniqueChunkData = unique(rawchunkData, 'rows');
    bias_ch2=mean(uniqueChunkData(:,3));%ch1 is timestamp, ch2 TMS, ch3 first muscle etc
    bias_ch3=mean(uniqueChunkData(:,4));%ch1 is timestamp, ch2 TMS, ch3 first muscle etc
    bias_ch4=mean(uniqueChunkData(:,5));%ch1 is timestamp, ch2 TMS, ch3 first muscle etc
    bias_ch5=mean(uniqueChunkData(:,6));%ch1 is timestamp, ch2 TMS, ch3 first muscle etc    
     
for ii = 1:nMEPs
    
    rawEpochData = [];
    
    triggerSent = 0;
    soundOccured=0;
    trialState = 1;
    exitTrial = false;
    flag=0;
    trialStartTime = GetSecs;   
    while exitTrial == false
        serialPortObj.TimerFcn = {'Rapid2_MaintainCommunication'};
        % Check the keyboard to see if a button has been pressed
        [keyIsDown,secs,keyCode] = KbCheck;
        
        % check for exit request
        if keyCode(escapeKey)
            exitTrial = true;
        end
        
        %  create variables from globals
        timeStamps = tS;
        rawData = rD;
        
        % send trigger if required
        if trialState == 3 && triggerSent == 0  ...
                && vbl - emgStartTime > triggerTime
            % send trigger via serial port here
            %fwrite(usb,marker); % write marker
            success = Rapid2_TriggerPulse(serialPortObj, 1)
            triggerSent = 1;
            mepCount = mepCount + 1;
            %clear old bytes from the buffer
            N = serialPortObj.BytesAvailable();
            while(N~=0)
                fread(serialPortObj,N);
                N = serialPortObj.BytesAvailable();
            end
        end
        
        
        
        
        %Caclulate the position of the feedback bar
        barYpos=screenYpixels-round(screenYpixels*(mean_MEP/Range));
        bar=[centreX-200 barYpos centreX+200 screenYpixels];
        
        %I include this just so a massive MEP doesnt crash the program if
        %it tries to plot it off the screen
        if barYpos<0
            barYpos=0;
        end
        
        if MEP_amp>mean_MEP
            barColour=[0 1 0];
        else
            barColour=[1 0 0];
        end
        
        %Caclulate the position of the rms background EMG bars
        
%         if FDI_rms<0.01
%             DotColourCh2_justbeforeMEP=[0 1 0];
%         else
%            DotColourCh2_justbeforeMEP =[1 0 0];
%         end
        
        
        
        % display settings for each trial state
        if trialState == 1  % rest
             rawRMSData = [rawRMSData; timeStamps rawData];
              uniqueRMSData = unique(rawRMSData, 'rows');
             RMS_chunk_ch2=uniqueRMSData(length(uniqueRMSData)-999:length(uniqueRMSData),3); %second channel is TMS. First is timestamp. Third is first EMG. Takes last 1000 samples for the chunk to take rms
             RMS_chunk_ch2=RMS_chunk_ch2-bias_ch2;
             rms_ch2=rms(RMS_chunk_ch2);
             RMS_chunk_ch3=uniqueRMSData(length(uniqueRMSData)-999:length(uniqueRMSData),4); %second channel is TMS. First is timestamp. Third is first EMG
             RMS_chunk_ch3=RMS_chunk_ch3-bias_ch3;
             rms_ch3=rms(RMS_chunk_ch3);
             RMS_chunk_ch4=uniqueRMSData(length(uniqueRMSData)-999:length(uniqueRMSData),5); %second channel is TMS. First is timestamp. Third is first EMG
             RMS_chunk_ch4=RMS_chunk_ch4-bias_ch4;
             rms_ch4=rms(RMS_chunk_ch4);
             RMS_chunk_ch5=uniqueRMSData(length(uniqueRMSData)-999:length(uniqueRMSData),6); %second channel is TMS. First is timestamp. Third is first EMG
             RMS_chunk_ch5=RMS_chunk_ch5-bias_ch5;
             rms_ch5=rms(RMS_chunk_ch5);
             
             if rms_ch2<backgnd_EMG_threshold
                  DotColourCh2=[0 1 0];
                  FDIrmsgood=1;
                 
             else
                 DotColourCh2=[1 0 0];
                 FDIrmsgood=0;
             end
             
             if rms_ch3<backgnd_EMG_threshold
                  DotColourCh3=[0 1 0];
                  ch3rmsgood=1;
             else
                 DotColourCh3=[1 0 0];
                 ch3rmsgood=0;
             end
             
             if rms_ch4<backgnd_EMG_threshold
                  DotColourCh4=[0 1 0];
                  ch4rmsgood=1;
             else
                 DotColourCh4=[1 0 0];
                 ch4rmsgood=0;
             end
             
             if rms_ch5<backgnd_EMG_threshold
                  DotColourCh5=[0 1 0];
                  ch5rmsgood=1;
             else
                 DotColourCh5=[1 0 0];
                 ch5rmsgood=0;
             end
             
             
             Screen('DrawDots', window, [centreX-1.5*dotSizePix centreY], dotSizePix, DotColourCh2, [], 2);  % EMG dot
             Screen('DrawDots', window, [centreX-0.5*dotSizePix centreY], dotSizePix, DotColourCh3, [], 2);  % EMG dot
             Screen('DrawDots', window, [centreX+0.5*dotSizePix centreY], dotSizePix, DotColourCh4, [], 2);  % EMG dot
             Screen('DrawDots', window, [centreX+1.5*dotSizePix centreY], dotSizePix, DotColourCh5, [], 2);  % EMG dot
            
        elseif trialState == 2 %fixation cross cue on screen
            rawRMSData = [rawRMSData; timeStamps rawData];
              uniqueRMSData = unique(rawRMSData, 'rows');
             RMS_chunk_ch2=uniqueRMSData(length(uniqueRMSData)-999:length(uniqueRMSData),3); %second channel is TMS. First is timestamp. Third is first EMG. Takes last 1000 samples for the chunk to take rms
             RMS_chunk_ch2=RMS_chunk_ch2-bias_ch2;
             rms_ch2=rms(RMS_chunk_ch2);
             RMS_chunk_ch3=uniqueRMSData(length(uniqueRMSData)-999:length(uniqueRMSData),4); %second channel is TMS. First is timestamp. Third is first EMG
             RMS_chunk_ch3=RMS_chunk_ch3-bias_ch3;
             rms_ch3=rms(RMS_chunk_ch3);
             RMS_chunk_ch4=uniqueRMSData(length(uniqueRMSData)-999:length(uniqueRMSData),5); %second channel is TMS. First is timestamp. Third is first EMG
             RMS_chunk_ch4=RMS_chunk_ch4-bias_ch4;
             rms_ch4=rms(RMS_chunk_ch4);
             RMS_chunk_ch5=uniqueRMSData(length(uniqueRMSData)-999:length(uniqueRMSData),6); %second channel is TMS. First is timestamp. Third is first EMG
             RMS_chunk_ch5=RMS_chunk_ch5-bias_ch5;
             rms_ch5=rms(RMS_chunk_ch5);
             
             if rms_ch2<backgnd_EMG_threshold
                  FDIrmsgood=1;
             else
                 DotColourCh2=[1 0 0];
                 FDIrmsgood=0;
             end
             
             if rms_ch3<backgnd_EMG_threshold
                  DotColourCh3=[0 1 0];
                  ch3rmsgood=1;
             else
                 DotColourCh3=[1 0 0];
                 ch3rmsgood=0;
             end
             
             if rms_ch4<backgnd_EMG_threshold
                  DotColourCh4=[0 1 0];
                  ch4rmsgood=1;
             else
                 DotColourCh4=[1 0 0];
                 ch4rmsgood=0;
             end
             
             if rms_ch5<backgnd_EMG_threshold
                  DotColourCh5=[0 1 0];
                  ch5rmsgood=1;
             else
                 DotColourCh5=[1 0 0];
                 ch5rmsgood=0;
             end
             
            Screen('DrawTexture', window, fixationTexture, [],[fixation_targetX-(fixationX/2) fixation_targetY-(fixationY/2) fixation_targetX+(fixationX/2) fixation_targetY+(fixationY/2)]);
        elseif trialState == 3  % emg
            rawEpochData = [rawEpochData; timeStamps rawData]; %aquiring the data to use for MEP measurement
            Screen('DrawTexture', window, fixationTexture, [],[fixation_targetX-(fixationX/2) fixation_targetY-(fixationY/2) fixation_targetX+(fixationX/2) fixation_targetY+(fixationY/2)]);
        elseif trialState == 4  % show feedback
            Screen('FillRect', window, barColour, bar);
            Screen('FillRect', window, lineColour, line);
%             Screen('FillRect', window, rmslineColour, rmsline);
%             Screen('FillRect', window, rms_bar_colour, FDIbar);
%             Screen('FillRect', window, rms_bar_colour, Muscle2bar);
%             Screen('FillRect', window, rms_bar_colour, Muscle3bar);
%             Screen('FillRect', window, rms_bar_colour, Muscle4bar);
%              Screen('DrawDots', window, [1000 960], dotSizePix, DotColourCh2_justbeforeMEP, [], 2);  % EMG dot
%              Screen('DrawDots', window, [1080 960], dotSizePix, DotColourCh2_justbeforeMEP, [], 2);  % EMG dot
%              Screen('DrawDots', window, [1160 960], dotSizePix, DotColourCh2_justbeforeMEP, [], 2);  % EMG dot
%              Screen('DrawDots', window, [1240 960], dotSizePix, DotColourCh2_justbeforeMEP, [], 2);  % EMG dot
           
             if MEP_amp>mean_MEP
                Screen('DrawTexture', window, tickTexture, [],[targetX-(tickX/2) targetY-(tickY/2) targetX+(tickX/2) targetY+(tickY/2)]);
                Screen('DrawTexture', window, dollarTexture, [],[dollar_targetX-(dollarX/2) dollar_targetY-(dollarY/2) dollar_targetX+(dollarX/2) dollar_targetY+(dollarY/2)]);
            DrawFormattedText(window,[num2str(nSuccess) '/' num2str(mepCount)], centreX+350, 840, white);
             else
                Screen('DrawTexture', window, crossTexture, [],[targetX-(crossX/2) targetY-(crossY/2) targetX+(crossX/2) targetY+(crossY/2)]);
            end
        end
        
        if MEP_amp>mean_MEP && flag==0 && trialState==4
            nSuccess=nSuccess+1;
            flag=1;
        end
       
        % Flip to the screen
        vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        
        % write results file
        writeData = [timeStamps, repmat(ii,[size(timeStamps),1]) ...
            rawData, repmat(trialState,[size(timeStamps),1])]';
        fprintf(dataFilePointer, [repmat('%i ', 1, size(writeData,1)-1) '%i\n'], writeData);
        
        % determine trial state
        if trialState == 1 && vbl - trialStartTime > gapDurationArray(ii,1) && FDIrmsgood==1 && ch3rmsgood==1 && ch4rmsgood==1 && ch5rmsgood==1 ;
            trialState = 2;
            fixStartTime = GetSecs;
            elseif trialState == 2 && vbl - fixStartTime > fixgapDurationArray(ii,1) && FDIrmsgood==1  && ch3rmsgood==1 && ch4rmsgood==1 && ch5rmsgood==1 ;
            success = Rapid2_ArmStimulator(serialPortObj);
            Rapid2_Delay(1000, serialPortObj);
            emgStartTime = GetSecs;
            trialState = 3;
        elseif trialState  == 3 && vbl - emgStartTime > emgDuration 
            trialState = 4;
            feedbackStartTime = GetSecs;
            % process emg data for feedback
            uniqueEpochData = unique(rawEpochData, 'rows');
            feedbackData = uniqueEpochData(:,2:end);
            save([direc filesep 'data' filesep subCode filesep SesNum filesep 'exp_UP' num2str(blockNumber) '-mep' (num2str(ii)) '.mat'], 'feedbackData','rawRMSData' )
            MEP_segment=feedbackData(1000:1600,2);%first EMG channel,ch2 on NI board
           
            %MEP_amp=peak2peak(MEP_segment);
            %For the control condition, the MEP amp is instead taken from
            %the created vector of amplitudes that were created earlier
            MEP_amp=control_MEP_vector(ii);
            
            %calculate the offset or bias, before rms measurement
            bias_ch2=mean(feedbackData(100:500,2));
            bias_ch3=mean(feedbackData(100:500,3));
            bias_ch4=mean(feedbackData(100:500,4));
            bias_ch5=mean(feedbackData(100:500,5));
            %remove bias and calculate rms background EMG
            biasCorrectedCh2=feedbackData(:,2)-bias_ch2;
            biasCorrectedCh3=feedbackData(:,3)-bias_ch3;
            biasCorrectedCh4=feedbackData(:,4)-bias_ch4;
            biasCorrectedCh5=feedbackData(:,5)-bias_ch5;
            
            
            %Provide auditory feedback
            if MEP_amp>mean_MEP && soundOccured ==0 && trialState == 4;
                PsychPortAudio('FillBuffer', pahandle, winSoundData);
                PsychPortAudio('Start', pahandle, 1, 0, 0);
                soundOccured=1;
            elseif MEP_amp<mean_MEP && soundOccured ==0;
                PsychPortAudio('FillBuffer', pahandle, loseSoundData);
                PsychPortAudio('Start', pahandle, 1, 0, 0);
                soundOccured=1;
            end
            
            %%% End while loop when time is up
        elseif trialState == 4 && vbl - feedbackStartTime > feedbackDuration
            exitTrial = true;
        end
        
        
    end
    
    % check for exit request
    if keyCode(escapeKey)
        break
    end
    
end

percentage_success=(nSuccess/nMEPs)*100
save([direc filesep 'data' filesep subCode filesep SesNum filesep 'UP_b' num2str(blockNumber) '-success' (num2str(nSuccess)) '.mat'], 'nSuccess', 'percentage_success')%save the num of successful trials
            

    function rData(src,event)
        rD = event.Data;
        tS = event.TimeStamps;
    end
%Clear serial port
% fclose(usb);
% delete(usb)           
% clear usb
delete(instrfindall);
PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
fclose('all');
close all;
clear all;

end