function recruitment_curve
%This version is currently suitable for the wireless delsys only.
%Plug the TMS into channel 1 and the EMG channels following this
clear all
ListenChar(0); %to ensure that the keyboard is responsive
rmt= input('What is the resting motor threshold?');
% set experiment parameters
subNum= input('Subject number');
subCode = ['NF' num2str(subNum)];
SessNumber = input('Session number');
direc = 'C:\Users\big lab\Documents\MATLAB';
%% file handling
subDir = ['C:\Users\big lab\Documents\MATLAB\data' filesep subCode filesep 'Recruitment_curves'];
if ~exist(subDir,'dir')
    mkdir(subDir);
end

dataFilename = [direc filesep 'data' filesep subCode filesep 'Recruitment_curves' filesep 'RC_ses'  num2str(SessNumber) '.dat'];

%This is protection for overwriting files. Keep commented for programming
%purposes but uncomment during real experiment
if fopen(dataFilename,'rt') ~= -1
    fclose('all');
    error('computer says no: result data file already exists!');
else
dataFilePointer = fopen(dataFilename,'wt');
end


nChans = 5;  %number of EMG channels
nMEPs = 60;  % number of MEPs to collect
tic
ninetyPercent= round(rmt*0.9);
hundredPercent= rmt;
hundredTenPercent= round(rmt*1.1);
hundredTwentyPercent= round(rmt*1.2);
hundredThirtyPercent= round(rmt*1.3);
hundredFortyPercent= round(rmt*1.4);
hundredFiftyPercent= round(rmt*1.5);
hundredSixtyPercent= round(rmt*1.6);
hundredSeventyPercent= round(rmt*1.7);
hundredEightyPercent= round(rmt*1.8);
IntensityVector=zeros(60,1);
IntensityVector(1:6,1)=ninetyPercent;
IntensityVector(7:12,1)=hundredPercent;
IntensityVector(13:18,1)=hundredTenPercent;
IntensityVector(19:24,1)=hundredTwentyPercent;
IntensityVector(25:30,1)=hundredThirtyPercent;
IntensityVector(31:36,1)=hundredFortyPercent;
IntensityVector(37:42,1)=hundredFiftyPercent;
IntensityVector(43:48,1)=hundredSixtyPercent;
IntensityVector(49:54,1)=hundredSeventyPercent;
IntensityVector(55:60,1)=hundredEightyPercent;

RandIntensityVector=IntensityVector(randperm(length(IntensityVector)));

%initialise empty vectors to store MEP amplitudes into
ninetyPercentMEPvector=[];
hundredPercentMEPvector=[];
hundredTenPercentMEPvector=[];
hundredTwentyPercentMEPvector=[];
hundredThirtyPercentMEPvector=[];
hundredFortyPercentMEPvector=[];
hundredFiftyPercentMEPvector=[];
hundredSixtyPercentMEPvector=[];
hundredSeventyPercentMEPvector=[];
hundredEightyPercentMEPvector=[];

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
Priority(2);
% The avaliable keys to press
spaceKey = KbName('SPACE');
% Dummy calls to make sure functions are ready to go without delay
KbCheck;
[keyIsDown, secs, keyCode] = KbCheck;


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
delete(instrfindall);  % to clear any pre existing COM port activities
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
success = Rapid2_ArmStimulator(serialPortObj);

   
%% START EXPERIMENT
for ii = 1:nMEPs
    disp(['MEP number ' num2str(ii)]);
%%set power level
powerLevel=RandIntensityVector(ii,1);

if powerLevel>100;
    powerLevel=100;
end

success = Rapid2_SetPowerLevel(serialPortObj, powerLevel, 1);
if ~success
    display 'Error: Cannot set the power level';
    return
else
    % Display power level;
    %display(powerLevel);
end
% Introduce delay to allow the  stimulator to adjust to the new power level
    Rapid2_Delay(1000, serialPortObj);
    
   %Pause to impose random time jitter between successive TMS pulses
    pause(gapDurationArray(ii,1));
%%   
    
switch powerLevel;
    case ninetyPercent
        intensityString='_90';
        
    case hundredPercent
        intensityString='_100';
        
    case hundredTenPercent
        intensityString='_110';
        
    case hundredTwentyPercent
        intensityString='_120';
        
    case hundredThirtyPercent
        intensityString='_130';
        
    case hundredFortyPercent
        intensityString='_140';
        
    case hundredFiftyPercent
        intensityString='_150';
        
    case hundredSixtyPercent
        intensityString='_160';
        
    case hundredSeventyPercent
        intensityString='_170';
        
    case hundredEightyPercent
        intensityString='_180';
end





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
    save(['C:\Users\big lab\Documents\MATLAB\data' filesep subCode filesep 'Recruitment_curves' filesep 'RC_ses'  num2str(SessNumber) '_int' intensityString '_' (num2str(ii)) '.mat'], 'feedbackData');
    
    TMS=uniqueEpochData(:,2);
    EMG1=uniqueEpochData(:,3);
    EMG2=uniqueEpochData(:,4);
    
    %Calculate MEP amplitude
    MEP_amp_EMG1=(max(EMG1(1245:1460)))-(min(EMG1(1245:1460)));
    MEP_amp_EMG2=(max(EMG2(1245:1460)))-(min(EMG2(1245:1460)));
    MEPvector1= [MEPvector1; MEP_amp_EMG1];
    MEPvector2= [MEPvector2; MEP_amp_EMG2];
    MeanMEP1=mean(MEPvector1); MeanMEP2=mean(MEPvector2);
    
    %Put the MEP amplitude measurment in the appropriate vector
    %corresponding to the TMS intensity
  switch powerLevel;
    case ninetyPercent
        ninetyPercentMEPvector=[ninetyPercentMEPvector; MEP_amp_EMG1];
        
    case hundredPercent
        hundredPercentMEPvector=[hundredPercentMEPvector; MEP_amp_EMG1];
        
    case hundredTenPercent
        hundredTenPercentMEPvector=[hundredTenPercentMEPvector; MEP_amp_EMG1];
        
    case hundredTwentyPercent
        hundredTwentyPercentMEPvector=[hundredTwentyPercentMEPvector; MEP_amp_EMG1];
        
    case hundredThirtyPercent
        hundredThirtyPercentMEPvector=[hundredThirtyPercentMEPvector; MEP_amp_EMG1];
        
    case hundredFortyPercent
        hundredFortyPercentMEPvector=[hundredFortyPercentMEPvector; MEP_amp_EMG1];
        
    case hundredFiftyPercent
        hundredFiftyPercentMEPvector=[hundredFiftyPercentMEPvector; MEP_amp_EMG1];
        
    case hundredSixtyPercent
        hundredSixtyPercentMEPvector=[hundredSixtyPercentMEPvector; MEP_amp_EMG1];
        
      case hundredSeventyPercent
        hundredSeventyPercentMEPvector=[hundredSeventyPercentMEPvector; MEP_amp_EMG1];
        
      case hundredEightyPercent
        hundredEightyPercentMEPvector=[hundredEightyPercentMEPvector; MEP_amp_EMG1];
end


    
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
    
    

    delete(EMGtext1);
    delete(MEPtext1);
    delete(EMGtext2);
    delete(MEPtext2);
    delete(MeanMEPtext1);
    delete(MeanMEPtext2);
    delete(plot1);
    delete(plot2);
end

meanNinety=mean(ninetyPercentMEPvector);
meanHundred=mean(hundredPercentMEPvector);
meanHundredTen=mean(hundredTenPercentMEPvector);
meanHundredTwenty=mean(hundredTwentyPercentMEPvector);
meanHundredThirty=mean(hundredThirtyPercentMEPvector);
meanHundredForty=mean(hundredFortyPercentMEPvector);
meanHundredFifty=mean(hundredFiftyPercentMEPvector);
meanHundredSixty=mean(hundredSixtyPercentMEPvector);
meanHundredSeventy=mean(hundredSeventyPercentMEPvector);
meanHundredEighty=mean(hundredEightyPercentMEPvector);
all_means=[meanNinety meanHundred meanHundredTen meanHundredTwenty meanHundredThirty meanHundredForty meanHundredFifty meanHundredSixty meanHundredSeventy meanHundredEighty];
hold on
plot(all_means);
title('Recruitment curve. Find which intensity corresponds to median MEP')
xlabel('TMS intensity level') % x-axis label
ylabel('MEP amplitude') % y-axis label
middle_MEP_size=median(all_means);
disp(['The median MEP amplitude was ' num2str(middle_MEP_size)]);
save(['C:\Users\big lab\Documents\MATLAB\data' filesep subCode filesep 'Recruitment_curves' filesep 'WholeRC_ses'  num2str(SessNumber) '.mat'], 'all_means');
    

    function rData(src,event)
        rD = event.Data;
        tS = event.TimeStamps;
    end
%Clear serial port
delete(instrfindall);

pause(10);
fclose('all');
close all;
%clear all;
toc
end