% REACTION TIME ANALYSIS AND ACCURACY ANALYSIS - detection
% Maria Ribeiro Feb 2016 updated June 2016 updated September 2016
% imperative times - reaction times go trials
% anticipatory_times - responses to cue (warning stimulus)
% slow responses with error warning imperative_times > 0.7 
% blank_times - trials with no imperative stimulus
% blank_keys - incorrectly pressed keys in blank trials

% clear
% Subjects_D_RT=zeros(14, 1);
participants=86; %[82 83 84 85];%[79 80 81]%77%78;%65:76;%62:64; %59:61; %[49:58];%[46 47 48]; %35;%[20 22]; %17;%31; %15;  %6; %[39 41:45]%[4 7:14 16 21 23 25:28 32:34 36:39 41:45];
task={'D1', 'D2', 'G1', 'G2'};
b=0;

for s=participants
b=b+1;
clearvars -except Subjects_D_RT task s Correlation_RT_Timing_Detection b

directory=strcat('M:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(s));

    for t=1:2; % 2 runs of detection task

    % load reaction times
    filename=strcat(directory, '\Behavior\', task(t), '\imperative_times.mat');
    load(filename{1});   
        if size(imperative_times, 2)>1;
            imperative_times(:,2:end)=[];
        else
        end
    change_var_name=strcat(task(t), '_RT=imperative_times;');
    eval(change_var_name{1});
    clear imperative_times;

    % load stimuli order
    filename=strcat(directory, '\Behavior\', task(t), '\stimuli_order.mat');
    load(filename{1});
    change_var_name=strcat(task(t), '_stimuli_order=stimuli_order;');
    eval(change_var_name{1});
    clear stimuli_order;
    
    % load impulsive responses to cue (warning stimulus)
    filename=strcat(directory, '\Behavior\', task(t), '\anticipatory_keys.mat');
    load(filename{1});
    change_var_name=strcat(task(t), '_anticipatory_keys=anticipatory_keys;');
    eval(change_var_name{1});
    clear anticipatory_keys;

    % load impulsive responses during blank trials
    filename=strcat(directory, '\Behavior\', task(t), '\blank_keys.mat');
    load(filename{1});
    change_var_name=strcat(task(t), '_blank_keys=blank_keys;');
    eval(change_var_name{1});
    clear blank_keys;
    
    %load Warning2ImperativeTime - time between warning and imperative
    %stimuli
    filename=strcat(directory, '\Behavior\', task(t), '\Warning2ImperativeTime.mat');
    load(filename{1});
    change_var_name=strcat(task(t), '_Warning2ImperativeTime=Warning2ImperativeTime;');
    eval(change_var_name{1});
    clear Warning2ImperativeTime;
    end

    % detection task
    %% checking for go trials
D1_GoTrials=find(D1_stimuli_order==1);
D2_GoTrials=find(D2_stimuli_order==1);

%% checking for trials where the participant incorrectly responded to the
% warning stimulus - anticipatory_responses==0;
if s==8;
    D1_responses2cue=cellfun('isempty',D1_anticipatory_keys);
    D1_responses2cue_trials=find(D1_responses2cue==0);

    D2_responses2cue=cellfun('isempty',D2_anticipatory_keys)';
    D2_responses2cue_trials=find(D2_responses2cue==0); 
    
else
    D1_responses2cue=cellfun('isempty',D1_anticipatory_keys);
    D1_responses2cue_trials=find(D1_responses2cue==0);

    D2_responses2cue=cellfun('isempty',D2_anticipatory_keys);
    D2_responses2cue_trials=find(D2_responses2cue==0); 
end
    
%% checking for trials with slow responses
D1_slow_responses=find(D1_RT>0.7);
D2_slow_responses=find(D2_RT>0.7);

%% checking for incorrect responses in blank trials
D1_responses2blanks_trials=[];
for n=1:length(D1_blank_keys);
    if ~isempty(D1_blank_keys{n});
        D1_responses2blanks_trials=[D1_responses2blanks_trials; n];
    end
end

D2_responses2blanks_trials=[];
for n=1:length(D2_blank_keys);
    if ~isempty(D2_blank_keys{n});
        D2_responses2blanks_trials=[D2_responses2blanks_trials; n];
    end
end
%% checking for missed responses
D1_MissedResponses=[];
for t=1:size(D1_RT,1);
    if D1_RT(t)==0 && D1_stimuli_order(t)==1;
D1_MissedResponses=[D1_MissedResponses; t];
    end
end
D2_MissedResponses=[];
for t=1:size(D2_RT,1);
    if D2_RT(t)==0 && D2_stimuli_order(t)==1;
D2_MissedResponses=[D2_MissedResponses; t];
    end
end

%% excluding incorrect trials + the trials immediately after
% does not exclude slow trials but excludes trials after slow trial given
% the error feedback that is presented and possibly might disrupt the next
% trial response
D1_incorrect_trials=cat(1, D1_responses2cue_trials,D1_responses2cue_trials+1, D1_slow_responses+1,...
     D1_MissedResponses, D1_MissedResponses+1, D1_responses2blanks_trials, D1_responses2blanks_trials+1);
D2_incorrect_trials=cat(1, D2_responses2cue_trials,D2_responses2cue_trials+1, D2_slow_responses+1,...
     D2_MissedResponses, D2_MissedResponses+1, D2_responses2blanks_trials, D2_responses2blanks_trials+1);

D_RT=cat(1, D1_RT(setdiff(D1_GoTrials, D1_incorrect_trials)), D2_RT(setdiff(D2_GoTrials, D2_incorrect_trials)));

% including in variable the warning to imperative time to check its effect
% on reaction time
D_RT(:,2)=1.5+cat(1, D1_Warning2ImperativeTime(setdiff(D1_GoTrials, D1_incorrect_trials))', D2_Warning2ImperativeTime(setdiff(D2_GoTrials, D2_incorrect_trials))');


%% saving data
mkdir(strcat(directory,'\BehavioralAnalysis'));
cd(strcat(directory,'\BehavioralAnalysis'));

% saving data
save D_RT_WithSlow D_RT; % reaction time excluding errors and trials after errors!
save D1_MissedResponses D1_MissedResponses;
save D2_MissedResponses D2_MissedResponses;
save D1_responses2cue_trials D1_responses2cue_trials
save D2_responses2cue_trials D2_responses2cue_trials
save D1_slow_responses D1_slow_responses;
save D2_slow_responses D2_slow_responses;
save D1_responses2blanks_trials D1_responses2blanks_trials;
save D2_responses2blanks_trials D2_responses2blanks_trials;

% % all subjects data
% D_RT_median=median(D_RT);
% Subjects_D_RT(s-2,1)=D_RT_median;

% plot
figure; plot(D_RT(:,1), D_RT(:,2),'d');
filename=strcat('AB', num2str(s));
title(filename,'FontSize',16,'FontWeight','bold');
xlabel('Reaction time (s)','FontSize',16,'FontWeight','bold');
ylabel('Warning to imperative time (s)','FontSize',16,'FontWeight','bold')

[R, P]=corrcoef(D_RT(:,1), D_RT(:,2));
% [RHO,PVAL] = corr(D_RT(:,1), D_RT(:,2), 'type', 'Spearman')

Correlation_RT_Timing_Detection(b)=R(1,2);

end

% cd('M:\\ProjectAgingNeuromodulation\AuditoryResearch\BehavioralAnalysis')
% save Correlation_RT_Timing_Detection Correlation_RT_Timing_Detection

[h p]=ttest(Correlation_RT_Timing_Detection)



