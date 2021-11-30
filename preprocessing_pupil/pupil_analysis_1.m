% 14 June 2016 - Maria Ribeiro
% Pupil dilation analysis
clear
participant='AB7';
task={'W1', 'D1', 'D2', 'G1', 'G2'};
task_initial={'W1', 'D1', 'D2', 'G1', 'G2'};

directory=strcat('L:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\ET\');
cd(directory);

%% transforming pupil data into .mat format
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for t=1%1:5;
filename_load=strcat(participant, '_', task_initial(t), ' Samples.txt');
filename_save=strcat(participant, '_', task_initial(t));
ET = parsesmi(filename_load{1} ,filename_save{1});
end

%% create variable with time and triggers for eye tracking data based on data from messages variable
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
% 12 - change in fixation cross for blinking

for t=1%1:5;
    clearvars -except t participant task task_initial directory
filename=strcat(participant, '_', task_initial(t));
load(filename{1})

event=[];
n=0;
messages_nowhite=strtrim(messages);

for i=1:length(messages_nowhite);    
    j=find(messages_nowhite{1,i}=='#');
            if strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: StartREC_detection');% warning only and detection task              
                n=n+1;
                event(n,2)=9;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
            elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: StarREC_gng');% gng task
                n=n+1;
                event(n,2)=9;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
            elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: inicrecordingposition');
                n=n+1;
                event(n,2)=8;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
            elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: EndRec');
                n=n+1;
                event(n,2)=10;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
            elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: cue');
                n=n+1;
                event(n,2)=1;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
            elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: target');
                n=n+1;
                event(n,2)=2;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
            elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: response');
                n=n+1;
                event(n,2)=5;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
             elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: SlowWarning');
                n=n+1;
                event(n,2)=3; 
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
             elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: nogo');
                n=n+1;
                event(n,2)=4;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
             elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: nogoresponse');
                n=n+1;
                event(n,2)=6;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
             elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: ImpulsiveResponseWarning');
                n=n+1;
                event(n,2)=7;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
             elseif strcmp(strtrim(messages_nowhite{1,i}((j+1):length(messages_nowhite{1,i}))),'Message: FixationChange');
                n=n+1;
                event(n,2)=12;
                event(n,1)=str2num(messages_nowhite{1,i}(1:j-8));
            end
end

colheader{1, length(colheader)}='Trigger';

% create column with events in data matrix
data(:,12)=0;
for e=1:size(event, 1);
    clear eventtime
    eventtime=find(data(:,1)<=event(e,1));
    data(eventtime(length(eventtime)), 12)=event(e,2);
end

% delete data before beginning of experiment
% % extract task
IniRec=find(data(:,12)==9);
% WarningData2=data(IniRec(1):end, :);
delete_trials=[1:IniRec(end)-1]; % only use the last task recorded
data(delete_trials,:)=[];
% extract task from event
IniRecEvent=find(event(:,2)==9);
% WarningEvents=event(IniRecEvent(end):end, :);
delete_events=[1:IniRecEvent(end)-1];
event(delete_events,:)=[];
% 
% delete trials after the end of task
EndRec=find(data(:,12)==10);
delete_end_trials=[EndRec(end)+1:size(data, 1)];
data(delete_end_trials,:)=[];
% extract task from event
EndRecEvent=find(event(:,2)==10); % changed column to 2! 26 october 2016!
% WarningEvents=event(IniRecEvent(end):end, :);
delete_end_events=[EndRecEvent(end)+1:size(event, 1)];
event(delete_end_events, :)=[];

filename_save=strcat(participant, '_', task(t), 'Events.mat');
% filename_save=strcat(participant, '_', 'D2Events.mat');
save(filename_save{1}, 'colheader', 'comments', 'data', 'messages', 'event');
end

