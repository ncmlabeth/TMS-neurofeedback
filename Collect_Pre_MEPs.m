
function Collect_Pre_MEPs
clear all
%This version is currently suitable for the wireless delsys only.
%Plug the TMS into channel 1 and the EMG channels following this
ListenChar(0); %to ensure that the keyboard is responsive
rmt= input('Desired Intensity');
% set experiment parameters
subNum= input('Subject number');
subCode = ['NF' num2str(subNum)];
SessionType= input('Is this an UP or DOWN session? Answer UP or DN    ', 's');
FirstOrSecondDay = input(['Is this the first or second session of the ' SessionType ' condition? Type 1 or 2      ']);
direc = 'C:\Users\big lab\Documents\MATLAB';
nChans = 5;  % number of EMG channels
nMEPs = 20;  % number of MEPs to collect

delete(instrfindall);
% if using matlab and psychtoolbox is installed, set to highest priority
% (9 for mac; 2 for windows)
% Priority(2);
% % initialise virtual com (serial) port
% usb = serial('COM3'); % virtual com port address will look something like this on windows
% % set necessary parameters; any program you use (eg presentation) should
% % allow you to specifiy these options
% set(usb,'BaudRate',115200);
% set(usb,'DataBits',8);
% set(usb,'Parity','none');
% set(usb,'StopBits',1);
% marker = 1; % the marker can be any value between 1 and 254
% % open serial port
% fopen(usb);
% % first write is always slow, so perform a dummy write before experiment
% % starts
% fwrite(usb,1);
%Pause to allow experimenter to get in postion
pause(10);
% Seed random number generator
rng('shuffle');
emgDuration = 1;  % collect EMG data for this duration
triggerTime = 0.5;  % time after trial start to send TMS trigger
minGap = 4; maxGap = 6;   % set min and max gap times
gapDurationArray = (rand(nMEPs,1)*maxGap) + minGap;  % array of gap duration between trials
MEPvector1=[];  MEPvector2=[];

%% Call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
% Give Matlab high priority
%Priority(2);
% The avaliable keys to press
spaceKey = KbName('SPACE');
% Dummy calls to make sure functions are ready to go without delay
KbCheck;
[keyIsDown, secs, keyCode] = KbCheck;

%% file handling
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

%dataFilename = ['data' filesep subCode filesep 'pre_' SessionType  num2str(FirstOrSecondDay) '-m' num2str(nMEPs) '-c' num2str(nChans) '.dat'];
dataFilename = [direc filesep 'data' filesep subCode filesep SesNum filesep 'pre_' SessionType  num2str(FirstOrSecondDay) '-m' num2str(nMEPs) '-c' num2str(nChans) '.dat'];

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

%Define some timings in frames for detecting MEP onset
triggerTimeFrames=triggerTime * s.Rate; %as 1000 frames is 0.5s with a 2000hz samplefreq
DelsysDelayFrames=0.048 * s.Rate;
MEP_latencyFrames=0.015 * s.Rate;
MEP_duration=0.045 * s.Rate;
%MEP_onsetFrames=triggerTimeFrames+DelsysDelayFrames+MEP_latencyFrames;
%MEP_offsetFrames=MEP_onsetFrames+MEP_duration;
%Define some timings in frames for background EMG rms calculation  
%backgroundEMG_start_relative_to_TMS=0.105 * s.Rate %start calculating 105ms before TMS
%backgroundEMG_end_relative_to_TMS=0.005 * s.Rate; %end calculating rms 5ms before TMS pulse
%backgroundEMG_startFrames=(triggerTimeFrames-backgroundEMG_start_relative_to_TMS)+DelsysDelayFrames;
%backgroundEMG_endFrames=(triggerTimeFrames-backgroundEMG_end_relative_to_TMS)+DelsysDelayFrames;
%% SET UP TMS AND SERIAL PORT COMMUNICATION
%delete(instrfindall);  % to clear any pre existing COM port activities
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
success = Rapid2_ArmStimulator(serialPortObj)
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
end
% Introduce delay to allow the  stimulator to adjust to the new power level
    Rapid2_Delay(1000, serialPortObj);
    pause(2);
%% START EXPERIMENT
for ii = 1:nMEPs
    
    rawEpochData = [];
    triggerSent = 0;
    exitTrial = false;
    trialStartTime = GetSecs;
    while exitTrial == false
        
        %  create variables from globals
        timeStamps = tS;
        rawData = rD;
        
        currentTime = GetSecs;
        rawEpochData = [rawEpochData; timeStamps rawData];
        % send trigger
        if triggerSent == 0 ...
                && currentTime - trialStartTime > triggerTime
            % send trigger via serial port here
%                fwrite(usb,marker); % write marker
            success = Rapid2_TriggerPulse(serialPortObj, 1)
            triggerSent = 1;
        end
        
        
        if  currentTime - trialStartTime > emgDuration
            exitTrial = true;
        end
        
        %This pause is key. Without this, the loop does not collect data
        %but runs anyway. Waitsecs doesnt work for this, and neither does a
        %while loop imposing timing using Getsecs > a criterion. Needs to
        %be pause. Also, enure the pause time is below the amount of time
        %needed to collect one bunch of samples, in this case 200 samples
        %in the listener so 100ms. Anything lower than this is ok. Larger
        %will mean there are gaps in the data acquisition
        pause(0.0001)
       
        %% Check the keyboard to see if a button has been pressed
        [keyIsDown,secs,keyCode] = KbCheck;
        % check for exit request
        if keyCode(spaceKey)
            MEPvector1=[];
            MEPvector1=[];
        end
    end
    
    %% process emg data for feedback
    uniqueEpochData = unique(rawEpochData, 'rows');
    feedbackData = uniqueEpochData(:,2:end);
    save([direc filesep 'data' filesep subCode filesep SesNum filesep 'pre_' SessionType num2str(FirstOrSecondDay) '-mep' (num2str(ii)) '.mat'], 'feedbackData');
    
    
    TMS=uniqueEpochData(:,2);
    EMG1=uniqueEpochData(:,3);
    EMG2=uniqueEpochData(:,4);
    
    %Calculate MEP amplitude
    MEP_amp_EMG1=(max(EMG1(1245:1460)))-(min(EMG1(1245:1460)));
    MEP_amp_EMG2=(max(EMG2(1245:1460)))-(min(EMG2(1245:1460)));
    MEPvector1= [MEPvector1; MEP_amp_EMG1];
    MEPvector2= [MEPvector2; MEP_amp_EMG2];
    MeanMEP1=mean(MEPvector1); MeanMEP2=mean(MEPvector2);
    
    %Calculate background EMG rms in 100 ms up until 5ms before TMS
    %calculate the offset or bias, before rms measurement
    bias_EMG1=mean(EMG1(100:500));
    bias_EMG2=mean(EMG2(100:500));
    %remove bias and calculate rms background EMG
    biasCorrectedEMG1=EMG1-bias_EMG1;
    biasCorrectedEMG2=EMG1-bias_EMG2;
    backgroundEMG1= rms(biasCorrectedEMG1(915:1115));
    backgroundEMG2= rms(biasCorrectedEMG2(915:1115));
    
    MEP_amp_string_EMG1= ['MEP ' num2str(MEP_amp_EMG1)];
    MeanMEP_string_EMG1= ['AVE ' num2str(MeanMEP1)];
    backgroundEMG_string= 'EMG';
    MEP_amp_string_EMG2= ['MEP ' num2str(MEP_amp_EMG2)];
    MeanMEP_string_EMG2= ['AVE ' num2str(MeanMEP2)];
    
    %Plotting
    subplot(2,1,1);
    ylim([-0.5 0.5]);
    xlim([1100 1450]);
    hold on
    %hline(0.05);
    %vline(triggerTimeFrames+DelsysDelayFrames);
    
    MEPtext1=text(1350,0.4, MEP_amp_string_EMG1, 'Color', 'red', 'FontSize', 34);
    MeanMEPtext1=text(1350,-0.2, MeanMEP_string_EMG1, 'Color', 'red', 'FontSize', 34);
    text(1110,0.2, backgroundEMG_string, 'Color', 'red', 'FontSize', 34);
    EMGtext1=text(1110,0.4, num2str(backgroundEMG1), 'Color', 'red', 'FontSize', 34);
    plot1=plot((biasCorrectedEMG1), 'LineWidth', 2);
    subplot(2,1,2);
    ylim([-0.5 0.5]);
    xlim([1100 1450]);
    hold on
    %hline(0.05);
    %vline(triggerTimeFrames+DelsysDelayFrames);
    MEPtext2=text(1350,0.4, MEP_amp_string_EMG2, 'Color', 'red', 'FontSize', 34);
    MeanMEPtext2=text(1350,-0.2, MeanMEP_string_EMG2, 'Color', 'red', 'FontSize', 34);
    text(1110,0.2, backgroundEMG_string, 'Color', 'red', 'FontSize', 34);
    EMGtext2=text(1110,0.4, num2str(backgroundEMG2), 'Color', 'red', 'FontSize', 34);
    plot2=plot((biasCorrectedEMG2), 'LineWidth', 2);
    hold off
    
    
    %Pause to impose random time jitter between successive TMS pulses
    pause(gapDurationArray(ii,1));
    delete(EMGtext1);
    delete(MEPtext1);
    delete(EMGtext2);
    delete(MEPtext2);
    delete(MeanMEPtext1);
    delete(MeanMEPtext2);
    delete(plot1);
    delete(plot2);
end

%Display the overall mean from the session, disregarding the first MEP
EndMeanMEP1=mean(MEPvector1(2:end));
EndSDMEP1=std(MEPvector1(2:end));
disp(['The mean MEP amplitude was ' num2str(EndMeanMEP1)]);
disp(['The standard deviation was ' num2str(EndSDMEP1)]);
    function rData(src,event)
        rD = event.Data;
        tS = event.TimeStamps;
    end
%Clear serial port
% fclose(usb);
% delete(usb)           
% clear usb
delete(instrfindall);
fclose('all');
close all;
%clear all;

end