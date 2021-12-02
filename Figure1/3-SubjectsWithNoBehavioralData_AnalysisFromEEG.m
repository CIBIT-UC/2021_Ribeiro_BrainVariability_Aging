% calculating behavioural variables for subjects from EEG data
% 24 Nov 2017 Maria Ribeiro

% trigger 11 - response after cue for both cue-only trials and correct
% trials

clear;
task_merged={'W1', 'D', 'G'};
% open eeglab 
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for p=[19 33 38];
    directory=strcat('E:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p),  '\EEG\');
    cd(directory);
    participant=strcat('AB', num2str(p));
% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
% load data detection - D1
filename=strcat(participant, '_D1_RefM1M2_RemEKG_EMG_Fps_filtered1_100Hz.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


% calculate behavioral variables  - detection task
clear events latencies
for e=1:size(EEG.event, 2);
    events(e, 1)=EEG.event(e).type;
    latencies(e, 1)=EEG.event(e).latency;
end

events([length(events)+1 length(events)+2  length(events)+3])=0; 

a=0;
for e=1:size(EEG.event, 2);
    if EEG.event(e).type==1;
        a=a+1;
        epoch_start(a, 1)=e;
    end
end

% behavioural variables
D1_RT=[]; D1_responses2cue_trials=[]; D1_slow_responses=[]; D1_responses2blanks_trials=[]; D1_MissedResponses=[]; 

for a=1:length(epoch_start);
    if isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 1]) || isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 10]); % correct trials
        D1_RT=[D1_RT; (latencies(epoch_start(a)+2)-latencies(epoch_start(a)+1))*0.001*2 (latencies(epoch_start(a)+1)-latencies(epoch_start(a)))*0.001*2];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 2]); % responses to cue
        D1_responses2cue_trials=[D1_responses2cue_trials; a]; 
    elseif  isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 3]); % slow trials
        D1_slow_responses=[D1_slow_responses; a];
        D1_RT=[D1_RT; (latencies(epoch_start(a)+3)-latencies(epoch_start(a)+2))*0.001*2 (latencies(epoch_start(a)+2)-latencies(epoch_start(a)+1))*0.001*2];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 1]) || isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 10]); % responses in cue-only trials
        D1_responses2blanks_trials=[D1_responses2blanks_trials; a];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2], 1),[2; 3]); % misses
        D1_MissedResponses=[D1_MissedResponses; a];
    end
end
        

%% load data detection - D2
filename=strcat(participant, '_D2_RefM1M2_RemEKG_EMG_Fps_filtered1_100Hz.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


% calculate behavioral variables  - detection task
clear events latencies epoch_start
for e=1:size(EEG.event, 2);
    events(e, 1)=EEG.event(e).type;
    latencies(e, 1)=EEG.event(e).latency;
end

events([length(events)+1 length(events)+2  length(events)+3])=0; 

a=0;
for e=1:size(EEG.event, 2);
    if EEG.event(e).type==1;
        a=a+1;
        epoch_start(a, 1)=e;
    end
end

% behavioural variables
D2_RT=[]; D2_responses2cue_trials=[]; D2_slow_responses=[]; D2_responses2blanks_trials=[]; D2_MissedResponses=[]; 

for a=1:length(epoch_start);
    if isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 1]) || isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 10]); ; % correct trials including trials where response was late
        D2_RT=[D2_RT; (latencies(epoch_start(a)+2)-latencies(epoch_start(a)+1))*0.001*2 (latencies(epoch_start(a)+1)-latencies(epoch_start(a)))*0.001*2];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 2]); % responses to cue
        D2_responses2cue_trials=[D2_responses2cue_trials; a]; 
    elseif  isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 3]); % slow trials
        D2_slow_responses=[D2_slow_responses; a];
         D2_RT=[D2_RT; (latencies(epoch_start(a)+3)-latencies(epoch_start(a)+2))*0.001*2 (latencies(epoch_start(a)+2)-latencies(epoch_start(a)+1))*0.001*2];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 1]) || isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 10]); % responses in cue-only trials
        D2_responses2blanks_trials=[D2_responses2blanks_trials; a];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2], 1),[2; 3]); % misses
        D2_MissedResponses=[D2_MissedResponses; a];
    end
end
      
D_RT=cat(1, D1_RT, D2_RT);
cd(strcat('E:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p),  '\BehavioralAnalysis\'));
save D_RT_WithSlow D_RT;
save D1_RT_WithSlow D1_RT;
save D2_RT_WithSlow D2_RT;
save D2_MissedResponses D2_MissedResponses;
save D2_responses2blanks_trials D2_responses2blanks_trials;
save D2_slow_responses D2_slow_responses;
save D2_responses2cue_trials D2_responses2cue_trials; 
save D1_MissedResponses D1_MissedResponses;
save D1_responses2blanks_trials D1_responses2blanks_trials;
save D1_slow_responses D1_slow_responses;
save D1_responses2cue_trials D1_responses2cue_trials;

%% gng task

    directory=strcat('E:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p),  '\EEG\');
    cd(directory);
    participant=strcat('AB', num2str(p));
% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
% load data gng - G1
filename=strcat(participant, '_G1_RefM1M2_RemEKG_EMG_Fps_filtered1_100Hz.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


% calculate behavioral variables  - detection task
clear events latencies
for e=1:size(EEG.event, 2);
    events(e, 1)=EEG.event(e).type;
    latencies(e, 1)=EEG.event(e).latency;
end


a=0;
for e=1:size(EEG.event, 2);
    if EEG.event(e).type==1;
        a=a+1;
        epoch_start(a, 1)=e;
    end
end

% behavioural variables
G1_RT=[]; G1_responses2cue_trials=[]; G1_slow_responses=[]; G1_responses2blanks_trials=[]; G1_MissedResponses=[]; 
G1_responses2nogo_trials=[];

for a=1:length(epoch_start);
    if isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 1]) || isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 10]); % correct trials
        G1_RT=[G1_RT; (latencies(epoch_start(a)+2)-latencies(epoch_start(a)+1))*0.001*2 (latencies(epoch_start(a)+1)-latencies(epoch_start(a)))*0.001*2];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 2]);% responses to cue
        G1_responses2cue_trials=[G1_responses2cue_trials; a]; 
    elseif  isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 3]); % slow trials
        G1_slow_responses=[G1_slow_responses; a];
        G1_RT=[G1_RT; (latencies(epoch_start(a)+3)-latencies(epoch_start(a)+2))*0.001*2 (latencies(epoch_start(a)+2)-latencies(epoch_start(a)+1))*0.001*2];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 1]) || isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 10]); % responses in cue-only trials
        G1_responses2blanks_trials=[G1_responses2blanks_trials; a];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2], 1),[2; 3]); % misses
        G1_MissedResponses=[G1_MissedResponses; a];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2], 1),[4; 6]);% responses in no-go trials
        G1_responses2nogo_trials=[G1_responses2nogo_trials; a]; 
    end
end
        

%% load data gng - G2
filename=strcat(participant, '_G2_RefM1M2_RemEKG_EMG_Fps_filtered1_100Hz.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


% calculate behavioral variables  - detection task
clear events latencies epoch_start
for e=1:size(EEG.event, 2);
    events(e, 1)=EEG.event(e).type;
    latencies(e, 1)=EEG.event(e).latency;
end

events([length(events)+1 length(events)+2  length(events)+3])=0; 

a=0;
for e=1:size(EEG.event, 2);
    if EEG.event(e).type==1;
        a=a+1;
        epoch_start(a, 1)=e;
    end
end

% behavioural variables
G2_RT=[]; G2_responses2cue_trials=[]; G2_slow_responses=[]; G2_responses2blanks_trials=[]; G2_MissedResponses=[]; 
G2_responses2nogo_trials=[];

for a=1:length(epoch_start);
    if isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 1]) || isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 10]); % correct trials including trials where response was late
        G2_RT=[G2_RT; (latencies(epoch_start(a)+2)-latencies(epoch_start(a)+1))*0.001*2 (latencies(epoch_start(a)+1)-latencies(epoch_start(a)))*0.001*2];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 2]); % responses to cue
        G2_responses2cue_trials=[G2_responses2cue_trials; a]; 
    elseif  isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[2; 5; 3]); % slow trials
        G2_slow_responses=[G2_slow_responses; a];
         G2_RT=[G2_RT; (latencies(epoch_start(a)+3)-latencies(epoch_start(a)+2))*0.001*2 (latencies(epoch_start(a)+2)-latencies(epoch_start(a)+1))*0.001*2];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 1]) || isequal(events([epoch_start(a)+1 epoch_start(a)+2 epoch_start(a)+3], 1),[5; 11; 10]); % responses in cue-only trials
        G2_responses2blanks_trials=[G2_responses2blanks_trials; a];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2], 1),[2; 3]); % misses
        G2_MissedResponses=[G2_MissedResponses; a];
    elseif isequal(events([epoch_start(a)+1 epoch_start(a)+2], 1),[4; 6]);% responses in no-go trials
        G2_responses2nogo_trials=[G2_responses2nogo_trials; a]; 
    end
end
      
G_RT=cat(1, G1_RT, G2_RT);
cd(strcat('E:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p),  '\BehavioralAnalysis\'));
save G_RT_WithSlow G_RT;
save G1_RT_WithSlow G1_RT;
save G2_RT_WithSlow G2_RT;
save G2_MissedResponses G2_MissedResponses;
save G2_responses2blanks_trials G2_responses2blanks_trials;
save G2_slow_responses G2_slow_responses;
save G2_responses2cue_trials G2_responses2cue_trials; 
save G1_MissedResponses G1_MissedResponses;
save G1_responses2blanks_trials G1_responses2blanks_trials;
save G1_slow_responses G1_slow_responses;
save G1_responses2cue_trials G1_responses2cue_trials; 
save G2_responses2nogo_trials G2_responses2nogo_trials;
save G1_responses2nogo_trials G1_responses2nogo_trials;
end