% Cued Auditory Task - detection condition
% Play different beeps prepare to respond and respond accordingly 
% Detection task - warning beep followed by imperative beep
% Go/no-go task - warning beep followed by target 80% of the times or no-go
% stimulus 20% of the time

% Clear the workspace
clear all;
close all;
sca;

% to make sure the random numbers are always different in different matlab
% sessions - warningtoimperative times and shuffle of trials
rand('twister',sum(100*clock))
randn('state',sum(100*clock))

% number of trials
Trials=60;
 
% variables to record
anticipatory_keys={};
anticipatory_times=[];
imperative_keys={}; 
imperative_times=[];
blank_keys=[];
blank_times=[];

% Number of channels and Frequency of the sound
nrchannels = 2;
freq = 192000;

% How many times do we wish to play the sound
repetitions = 1;
% Length of the beep
beepLengthSecs = 0.250;
% Length of the pause between warning beep and imperative beep
% Exponential random number using function exprnd - number of trial 100
% Warning2ImperativeTime=exprnd(0.25, 50, 1);
for i = 1:60
U = rand(1);
Warning2ImperativeTime(i) = -0.25*log(U);
end
% Length of the pause between imperative beep and warning beep
% Exponential random number using function exprnd - number of trial 50
for j = 1:60
U = rand(1);
Imperative2WarningTime(j) = -1*log(U);
end
% beepPauseTime = 1;
% Start immediately (0 = immediately)
startCue = 0;
% Should we wait for the device to really start (1 = yes)
% INFO: See help PsychPortAudio
waitForDeviceStart = 1;

% Make a beep for warning cue
WarningBeep = MakeBeep(1500, beepLengthSecs, freq)*0.2;

% Make Beep for imperative stimuli
TargetBeep= MakeBeep(1700, beepLengthSecs, freq)*0.1;

% Make Beep for slow responses
SlowResponsesBeep = MakeBeep(1000, beepLengthSecs, freq)*0.3;

% order of presentation of go and blank trials
% 1 - go stimuli; 3 - blank no imperative stimulus
stimuli_order=[];
target_flanking=[1 1];
blank_order=[1 1 1 3];

for g=1:10
    stimuli_order=[stimuli_order target_flanking Shuffle(blank_order)];
end

% to get out of Psychtoolbox in an error occurs
% try
KbQueueCreate; % for measuring key presses   
KbQueueStart; %start listening for key presses
KbName('UnifyKeyNames');
ListenChar(2);
    
% show grey background with fixation dot
    Screen('Preference','Verbosity',0)
    Screen('Preference','SkipSyncTests',1);
    Screen('Preference','VisualDebugLevel',0);
    
    HideCursor;
    
    screennumber=0;% screen to display stimuli
    [windowPtr,rect]=Screen('OpenWindow',screennumber); 
    monitorFlipInterval=Screen('GetFlipInterval', windowPtr);
    [X,Y] = RectCenter(rect);
%     resolucao=2*X;

    ivx.window = windowPtr;
    
    %% To send EEG triggers
    address=888;

%% Remote access to iViewX

    pnet('closeall');

    host='192.168.2.1';
    port=4444;
    ivx= iViewXInitDefaults(ivx.window, 5  , host, port); % no [] colocar o numero de pontos que quiser na calibração
    ivx.localport=4444;
    ivx;
 ivx.backgroundColour= [120 120 120];
    
    ivx.screenHSize =  2*X;
    ivx.screenVSize = 2*Y;

    % calibration
    [result, ivx]=iViewX('calibrate', ivx);

    
    %% blank the Screen
    Screen('FillRect',windowPtr,[120 120 120]);
    Screen(windowPtr, 'Flip');
           tcal=GetSecs;
    while GetSecs-tcal<1.5
        Screen('TextSize', windowPtr, 50);
        Screen('DrawText', windowPtr, 'Vamos começar', X-250, Y-45, [50 50 50]);
        vbl=Screen(windowPtr, 'Flip');
    end
    WaitSecs(1);

    % Fixation Cross
    FixCross = [X-1,Y-20,X+1,Y+20;X-20,Y-1,X+20,Y+1];
    Screen('FillRect', windowPtr, [140, 140, 140], FixCross');
    Screen('Flip', windowPtr);

  

        %% Set Recording
            % start recording eye tracking
            [result, ivx]=iViewX('startrecording', ivx);
            

%% wait for pupil to stabilize
WaitSecs(5);
            %trigger start recording
            %TRIGGERS EEG
            lptwrite(address, 9);
            WaitSecs(0.004);
            lptwrite(address, 0);
            iViewXComm('send', ivx, 'ET_AUX "StartREC_detection"');
 
WaitSecs(5); 
        
            
            %% Initialize Sounddriver
InitializePsychSound(1);

% Open Psych-Audio port, with the follow arguments
% (1) [] = default sound device
% (2) 1 = sound playback only
% (3) 1 = default level of latency
% (4) Requested frequency in samples per second
% (5) 2 = stereo output
pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels);

% Set the volume to half for this demo
PsychPortAudio('Volume', pahandle, 0.75);

% record experiment start time
start_time=GetSecs;



    for t=1:Trials; % loop through trials
        
KbQueueFlush(); %clear any pressed key
        
% Fill the audio playback buffer with the audio data, doubled for stereo
% presentation
% Warning stimulus - cue
PsychPortAudio('FillBuffer', pahandle, [WarningBeep; WarningBeep]);

% to count number of responses to stimulus
count=0;
key=[];

% Start audio playback
PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
cue_time(t)=GetSecs;

            %TRIGGERS EEG ET
            lptwrite(address, 1);
            WaitSecs(0.004);
            lptwrite(address, 0);
            iViewXComm('send', ivx, 'ET_AUX "cue"');%trigger eye tracker

% The beep will play for 0.250 second, so we have to wait for that length of
% time, PLUS the amount of time we want between beeps
% exponential distribution random generated numbers - minimum 1.5 sec 
while GetSecs-cue_time(t)<(1.5 + Warning2ImperativeTime(t));
     [keyIsDown, firstPress] = KbQueueCheck;
                        if keyIsDown==1

                            %TRIGGERS PARA A RESPOSTA
                            lptwrite(address, 5);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, strcat('ET_AUX "response"'));
                            
                            count=count+1;
                            pressedCodes=find(firstPress);
                            for j=1:size(pressedCodes,2)
                                key = KbName(pressedCodes(j));
                                secs=firstPress(pressedCodes(j))-cue_time(t);
                            end

                            anticipatory_keys{t, count}=key;
                            anticipatory_times(t, count)=secs;
                            
                            % play error sound
                             WaitSecs(0.2);
                            % Put sound in buffer   
    PsychPortAudio('FillBuffer', pahandle, [SlowResponsesBeep;SlowResponsesBeep]);
    % Start audio playback
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
                            %TRIGGERS EEG ET - error feedback
                            lptwrite(address, 11);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, 'ET_AUX "ErrorWarning"');%trigger eye tracker

                        end  
end

if stimuli_order(t)==1; % play go stimuli
% Fill the audio playback buffer with the audio data, doubled for stereo
% presentation
% Imperative stimulus - target
PsychPortAudio('FillBuffer', pahandle, [TargetBeep;TargetBeep]);

%count number of responses
count=0;
key=[];

KbQueueFlush(); %clear any pressed key

% Start audio playback
PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
target_time(t)=GetSecs;
                            %TRIGGERS EEG
                            lptwrite(address, 2);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, strcat('ET_AUX "target"'));%trigger eye tracker

% Measure response time and response key
while GetSecs-target_time(t)<1.2; %wait 5 secs before next stimulus
                        [keyIsDown, firstPress] = KbQueueCheck;
                        if keyIsDown==1

                            %TRIGGERS PARA A RESPOSTA
                            lptwrite(address, 5);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, strcat('ET_AUX "response"'));
                            
                            count=count+1;
                            pressedCodes=find(firstPress);
                            for j=1:size(pressedCodes,2)
                                key = KbName(pressedCodes(j));
                                secs=firstPress(pressedCodes(j))-target_time(t);
                            end
                            imperative_keys{t, count}=key;
                            imperative_times(t, count)=secs;
                        end
end

% play error sound for slow or missed responses
if size(imperative_times,1)<t; %missed responses
    % Put sound in buffer   
    PsychPortAudio('FillBuffer', pahandle, [SlowResponsesBeep;SlowResponsesBeep]);
 % Start audio playback
PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
                            %TRIGGERS EEG ET - error feedback
                            lptwrite(address, 3);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, 'ET_AUX "SlowWarning"');%trigger eye tracker
imperative_times(t, 1)=0;

elseif imperative_times(t, 1)>0.7; %slow responses
    % Put sound in buffer   
    PsychPortAudio('FillBuffer', pahandle, [SlowResponsesBeep;SlowResponsesBeep]);
    % Start audio playback
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
                            %TRIGGERS EEG ET - error feedback
                            lptwrite(address, 3);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, 'ET_AUX "SlowWarning"');%trigger eye tracker
end

elseif stimuli_order(t)==3; % do not play imperative stimulus - blank trial
target_time(t)=GetSecs;    
       %count number of responses
count=0;
key=[]; 
     while GetSecs-target_time(t)<1.2; %wait 1.2 secs for response
                        [keyIsDown, firstPress] = KbQueueCheck;
                        if keyIsDown==1

                            %TRIGGERS PARA A RESPOSTA
                            lptwrite(address, 5);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, strcat('ET_AUX "blankresponse"'));
                            
                            count=count+1;
                            pressedCodes=find(firstPress);
                            for j=1:size(pressedCodes,2)
                                key = KbName(pressedCodes(j));
                                secs=firstPress(pressedCodes(j))-target_time(t);
                            end
                            blank_keys{t, count}=key;
                            blank_times(t, count)=secs;
                            
                            % play error sound
                             WaitSecs(0.2);
                            % Put sound in buffer   
    PsychPortAudio('FillBuffer', pahandle, [SlowResponsesBeep;SlowResponsesBeep]);
    % Start audio playback
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
                            %TRIGGERS EEG ET - error feedback
                            lptwrite(address, 11);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, 'ET_AUX "ErrorWarning"');%trigger eye tracker

                        end
    end

        % play error sound for impulsive responses
        if size(blank_times,1)<t;
        blank_times(t, 1)=-1; % correct withhold of response

        elseif size(blank_times,1)==t; % subject incorrectly pressed response button - don't give feedback in these case
%         Put sound in buffer   
%         PsychPortAudio('FillBuffer', pahandle, [ImpulsiveResponsesBeep;ImpulsiveResponsesBeep]);
%         Start audio playback
%         PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
%         iViewXComm('send', ivx, 'ET_AUX "ImpulsiveResponseWarning"');%trigger eye tracker
        end    
        
end
temp_time=GetSecs;
% WaitSecs(4+Imperative2WarningTime(t));
     while GetSecs-temp_time<(4+Imperative2WarningTime(t)); %wait between trials
                        [keyIsDown, firstPress] = KbQueueCheck;
                        if keyIsDown==1

                            %TRIGGERS PARA A RESPOSTA
                            lptwrite(address, 5);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, strcat('ET_AUX "response"'));
                            
                            % play error sound
                             WaitSecs(0.2);
                            % Put sound in buffer   
    PsychPortAudio('FillBuffer', pahandle, [SlowResponsesBeep;SlowResponsesBeep]);
    % Start audio playback
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
                            %TRIGGERS EEG ET - error feedback
                            lptwrite(address, 11);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, 'ET_AUX "ErrorWarning"');%trigger eye tracker

                        end
    end

    end

   WaitSecs(10);
   
   % send end trigger
   
                            %TRIGGERS EEG
                            lptwrite(address, 10);
                            WaitSecs(0.004);
                            lptwrite(address, 0);
                            iViewXComm('send', ivx, 'ET_AUX "EndRec"');%trigger eye tracker
% Stop playback
PsychPortAudio('Stop', pahandle);
% Close the audio device
PsychPortAudio('Close', pahandle);

% catch
%     % This section is executed only in case an error happens in the
%     % experiment code implemented between try and catch...
%     ShowCursor;
%     Screen('CloseAll'); %or sca
%     ListenChar(0);
% %     Screen('Preference', 'VisualDebuglevel', olddebuglevel);
%     %output the error message
%     psychrethrow(psychlasterror);
% end
% end
     % stop recording eye tracking
     [result, ivx]=iViewX('stoprecording', ivx);
  
% create folder and save behavioral data     
mkdir(strcat(pwd, '\D2'));
cd(strcat(pwd, '\D2'));
 save imperative_keys imperative_keys
 save imperative_times imperative_times
 save anticipatory_keys anticipatory_keys
 save anticipatory_times anticipatory_times 
 save target_time target_time
 save blank_keys blank_keys
 save blank_times blank_times                           
 save Warning2ImperativeTime Warning2ImperativeTime
 save Imperative2WarningTime Imperative2WarningTime
save stimuli_order stimuli_order

    Screen('CloseAll');
    ShowCursor;
    
    
            % save recorded data - cannot save remote data so far!
%         user = 'User1';
%         description = 'Description1';
%         ovr = int32(1);
%         current_wd = pwd;
%         filename = fullfile(current_wd,['iViewXSDK_Matlab_Slideshow_Data_' user '.idf']);
%         ret_save = iView.iV_SaveData(filename, description, user, ovr);
%         
%         if(ret_save ~= 1)
%             disp('Recording data could not be saved')
%         end
%     
ListenChar(0);