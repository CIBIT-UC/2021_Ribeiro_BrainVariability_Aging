% behavioral analysis by eeg triggers - after rejection of periods with
% artifacts!
% Maria Ribeiro - May 2016
% Triggers 
% 9 - start recording StarRec
% 10 - end recording EndRec
% 8 - start recording inicrecordingposition
% 1 - cue
% 2 - target
% 3 - slow warning
% 4 - no-go stimulus
% 5 - response to target or cue
% 6 - impulsive response to no-go stimulus
% 7 - warning following impulsive response
% 11 - error warning after incorrect response to the cue or error responses between trials "ErrorWarning"
clear

younger=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];

participants=[younger older];
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
task_merged={'W1', 'D', 'G'};

for p=participants
participant=strcat('AB', num2str(p));
directory1=strcat('L:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\EEG\');
directory=strcat('L:\\ProjectAgingNeuromodulation\\AuditoryResearch\\AuditoryTask_EyeTracking_EEG_lab94\\', participant, '\\EEG\\');

cd(directory1);

% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

%load file in eeglab - detection task
t=2;
filename=strcat(participant, '_', task_merged{t}, '_Merged_epochs_baselined.set');
EEG = pop_loadset(filename);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

clear epoch_event;
epoch_event=zeros(size(EEG.event, 2), 3);

for e=1:size(EEG.event, 2)
    epoch_event(e, 1)=EEG.event(e).epoch;
    epoch_event(e, 2)=str2double(EEG.event(e).type);
    epoch_event(e, 3)=EEG.event(e).latency;
end

cue_only_correct_trial=[];
correct_trial=[];
response2cue=[];
multiple_responses=[];
misses=[];
response2cue=[];
detection_RT=[];
slow_responses=[];
slow_RT=[];
warning2imperativetime_detection=[];

for r=1:epoch_event(end, 1)
    events_per_epoch_index=[];
    events_per_epoch_index=find(epoch_event(:,1)==r);
    
    clear events_per_epoch latencies_per_epoch
    events_per_epoch=epoch_event(events_per_epoch_index, 2);
    latencies_per_epoch=epoch_event(events_per_epoch_index, 3);
    
    if length(events_per_epoch)==1
        cue_only_correct_trial=[cue_only_correct_trial; r];
    else% changed 17Jan2018 %if length(events_per_epoch)>=3;
        if events_per_epoch(2)==2 && events_per_epoch(3)==5 && length(events_per_epoch)==3
            correct_trial=[correct_trial; r];
            detection_RT=[detection_RT; r (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/500]; % corrected samplig rate 27April2018
            warning2imperativetime_detection=[warning2imperativetime_detection; r (latencies_per_epoch(2)-latencies_per_epoch(1))*1000/500]; % corrected samplig rate 27April2018
        elseif events_per_epoch(2)==2 && events_per_epoch(3)==5 && events_per_epoch(4)==3 % corrected 4 october 2016
            slow_responses=[slow_responses; r];
            slow_RT=[slow_RT; r (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/500];% corrected samplig rate 27April2018
        elseif events_per_epoch(2)==5
            response2cue=[response2cue; r];
        elseif length(events_per_epoch)>3 && events_per_epoch(2)==2 && events_per_epoch(3)==5 && events_per_epoch(4)==5
            multiple_responses=[multiple_responses; r]; 
        elseif events_per_epoch(2)==2 && events_per_epoch(3)==3 % changed 17Jan2018  
            misses=[misses; r];
        end
    end 
end


mkdir(strcat(pwd, '\detection_behavior'));
cd(strcat(pwd, '\detection_behavior'));

save cue_only_correct_trial cue_only_correct_trial
save response2cue response2cue
save correct_trial correct_trial
save detection_RT detection_RT
save multiple_responses multiple_responses
save misses misses
save slow_responses slow_responses
save slow_RT slow_RT
save warning2imperativetime_detection warning2imperativetime_detection

cd(directory1);

%% load file in eeglab - go/no-go task
t=3;
filename=strcat(participant, '_', task_merged{t}, '_Merged_epochs_baselined.set');
EEG = pop_loadset(filename);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

clear epoch_event;
epoch_event=zeros(size(EEG.event, 2), 3);

for e=1:size(EEG.event, 2)
    epoch_event(e, 1)=EEG.event(e).epoch;
    epoch_event(e, 2)= str2double(EEG.event(e).type);
    epoch_event(e, 3)=EEG.event(e).latency;
end

cue_only_correct_trial=[];
correct_trial=[];
response2cue=[];
multiple_responses=[];
misses=[];
response2cue=[];
correct_nogo=[];
response2nogo=[];
go_RT=[];
slow_responses=[];
slow_RT=[];
warning2imperativetime_gng=[];

for w=1:epoch_event(end, 1)
    events_per_epoch_index=[];
    events_per_epoch_index=find(epoch_event(:,1)==w);
    
    clear events_per_epoch latencies_per_epoch
    events_per_epoch=epoch_event(events_per_epoch_index, 2);
    latencies_per_epoch=epoch_event(events_per_epoch_index, 3);
    
    if length(events_per_epoch)==1
        cue_only_correct_trial=[cue_only_correct_trial; w];
    elseif length(events_per_epoch)==2 && events_per_epoch(2)==4% changed 17Jan2018   
        correct_nogo=[correct_nogo; w];   
    elseif length(events_per_epoch)>=3
        if events_per_epoch(2)==2 && events_per_epoch(3)==5 && length(events_per_epoch)==3
            correct_trial=[correct_trial; w];
            go_RT=[go_RT; w (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/500];% corrected samplig rate 27April2018
            warning2imperativetime_gng=[warning2imperativetime_gng; w (latencies_per_epoch(2)-latencies_per_epoch(1))*1000/500]; % corrected samplig rate 27April2018
        elseif events_per_epoch(2)==2 && events_per_epoch(3)==5 && events_per_epoch(4)==3 % corrected 4 october 2016
            slow_responses=[slow_responses; w];
            slow_RT=[slow_RT; w (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/500]; % corrected samplig rate 27April2018
        elseif events_per_epoch(2)==5
            response2cue=[response2cue; w];
        elseif length(events_per_epoch)>3 && events_per_epoch(2)==2 && events_per_epoch(3)==5 && events_per_epoch(4)==5
            multiple_responses=[multiple_responses; w];
        elseif events_per_epoch(2)==4 && events_per_epoch(3)==6 % changed 17Jan2018  
            response2nogo=[response2nogo; w];  
        elseif events_per_epoch(2)==2 && events_per_epoch(3)==3 % changed 17Jan2018  
            misses=[misses; w];
        end  
    end 
            
end

mkdir(strcat(pwd, '\GNG_behavior'));
cd(strcat(pwd, '\GNG_behavior'));

save cue_only_correct_trial cue_only_correct_trial
save response2cue response2cue
save correct_trial correct_trial
save multiple_responses multiple_responses
save misses misses
save correct_nogo correct_nogo
save response2nogo response2nogo
save slow_responses slow_responses
save slow_RT slow_RT
save go_RT go_RT
save warning2imperativetime_gng warning2imperativetime_gng

end