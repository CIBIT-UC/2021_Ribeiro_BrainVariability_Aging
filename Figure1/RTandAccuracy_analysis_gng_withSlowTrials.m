% REACTION TIME ANALYSIS AND ACCURACY ANALYSIS - gng
% Maria Ribeiro Feb 2016 updated June 2016
% imperative times - reaction times go trials
% anticipatory_times - responses to cue (warning stimulus)
% slow responses with error warning imperative_times > 0.7 
% blank_times - trials with no imperative stimulus
% blank_keys - incorrectly pressed keys in blank trials
% nogo_keys - incorrect responses to no-go stimulus

clear
younger=[4 9 10 13 15 16 25 26 28 31 34 36 42 44 45 46 50 51 53 54 56 59 62 66 68 72 74 76 78 80 81 82 84 85]; % 33 behaviour not saved
older=[7 8	11	12	14	17	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77  79 83   86]; % 19 behaviour not saved

Subjects_G_RT=zeros(11, 1);
participants=[younger older];
task={'D1', 'D2', 'G1', 'G2'};
b=0;
for s=participants;
b=b+1;
clearvars -except Subjects_G_RT task s participants b Correlation_RT_Timing_GNG

directory=strcat('M:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(s));

    for t=3:4; % 2 runs of gng task

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

    % load impulsive responses durign blank trials
    filename=strcat(directory, '\Behavior\', task(t), '\blank_keys.mat');
    load(filename{1});
    change_var_name=strcat(task(t), '_blank_keys=blank_keys;');
    eval(change_var_name{1});
    clear blank_keys;
    
    % load impulsive responses to no-go stimulus
    filename=strcat(directory, '\Behavior\', task(t), '\nogo_keys.mat');
    load(filename{1});
    change_var_name=strcat(task(t), '_nogo_keys=nogo_keys;');
    eval(change_var_name{1});
    clear nogo_keys;
    
    %load Warning2ImperativeTime - time between warning and imperative
    %stimuli
    filename=strcat(directory, '\Behavior\', task(t), '\Warning2ImperativeTime.mat');
    load(filename{1});
    change_var_name=strcat(task(t), '_Warning2ImperativeTime=Warning2ImperativeTime;');
    eval(change_var_name{1});
    clear Warning2ImperativeTime;
    end

    % gng task
    %checking for go trials
G1_GoTrials=find(G1_stimuli_order==1);
G2_GoTrials=find(G2_stimuli_order==1);

% checking for trials where the participant incorrectly responded to the
% warning stimulus - anticipatory_responses==0;
if s==8;
    G1_responses2cue=cellfun('isempty',G1_anticipatory_keys);
    G1_responses2cue_trials=find(G1_responses2cue==0);

    G2_responses2cue=cellfun('isempty',G2_anticipatory_keys)';
    G2_responses2cue_trials=find(G2_responses2cue==0); 
    
else
    G1_responses2cue=cellfun('isempty',G1_anticipatory_keys);
    G1_responses2cue_trials=find(G1_responses2cue==0);

    G2_responses2cue=cellfun('isempty',G2_anticipatory_keys);
    G2_responses2cue_trials=find(G2_responses2cue==0); 
end
    
%% checking for trials with slow responses
G1_slow_responses=find(G1_RT>0.7);
G2_slow_responses=find(G2_RT>0.7);

%% checking for incorrect responses in blank trials
G1_responses2blanks_trials=[];
for n=1:length(G1_blank_keys);
    if ~isempty(G1_blank_keys{n});
        G1_responses2blanks_trials=[G1_responses2blanks_trials; n];
    end
end

G2_responses2blanks_trials=[];
for n=1:length(G2_blank_keys);
    if ~isempty(G2_blank_keys{n});
        G2_responses2blanks_trials=[G2_responses2blanks_trials; n];
    end
end

%% checking for missed responses - it was wrong! corrected 18 August 2017
G1_MissedResponses=[];
for t=1:size(G1_RT,1);
    if G1_RT(t)==-1 && G1_stimuli_order(t)==1;
G1_MissedResponses=[G1_MissedResponses; t];
    end
end
G2_MissedResponses=[];
for t=1:size(G2_RT,1);
    if G2_RT(t)==-1 && G2_stimuli_order(t)==1;
G2_MissedResponses=[G2_MissedResponses; t];
    end
end

%% checking for incorrect responses to no-go stimulus
G1_responses2nogo_trials=[];
for n=1:length(G1_nogo_keys);
    if ~isempty(G1_nogo_keys{n});
        G1_responses2nogo_trials=[G1_responses2nogo_trials; n];
    end
end

G2_responses2nogo_trials=[];
for n=1:length(G2_nogo_keys);
    if ~isempty(G2_nogo_keys{n});
        G2_responses2nogo_trials=[G2_responses2nogo_trials; n];
    end
end

%% excluding incorrect trials + the trials immediately after
% does not exclude slow trials but excludes trials after slow trial given
% the error feedback that is presented and possibly might disrupt the next
% trial response
G1_incorrect_trials=cat(1, G1_responses2cue_trials,G1_responses2cue_trials+1, G1_slow_responses+1,...
     G1_MissedResponses, G1_MissedResponses+1, G1_responses2blanks_trials, G1_responses2blanks_trials+1, G1_responses2nogo_trials, G1_responses2nogo_trials+1);
G2_incorrect_trials=cat(1, G2_responses2cue_trials,G2_responses2cue_trials+1, G2_slow_responses+1,...
     G2_MissedResponses, G2_MissedResponses+1, G2_responses2blanks_trials, G2_responses2blanks_trials+1, G2_responses2nogo_trials, G2_responses2nogo_trials+1);

 % append reaction time from both runs only correct trials
G_RT(:,1)=cat(1, G1_RT(setdiff(G1_GoTrials, G1_incorrect_trials)), G2_RT(setdiff(G2_GoTrials, G2_incorrect_trials)));
% including in variable the warning to imperative time to check its effect
% on reaction time
G_RT(:,2)=1.5+cat(1, G1_Warning2ImperativeTime(setdiff(G1_GoTrials, G1_incorrect_trials))', G2_Warning2ImperativeTime(setdiff(G2_GoTrials, G2_incorrect_trials))');


%% saving data
% mkdir(strcat(directory,'\BehavioralAnalysis'));
cd(strcat(directory,'\BehavioralAnalysis'));

% saving data
save G_RT_WithSlow G_RT;
save G1_MissedResponses G1_MissedResponses;
save G2_MissedResponses G2_MissedResponses;
save G1_responses2cue_trials G1_responses2cue_trials
save G2_responses2cue_trials G2_responses2cue_trials
save G1_slow_responses G1_slow_responses;
save G2_slow_responses G2_slow_responses;
save G1_responses2blanks_trials G1_responses2blanks_trials;
save G2_responses2blanks_trials G2_responses2blanks_trials;
save G1_responses2nogo_trials G1_responses2nogo_trials 
save G2_responses2nogo_trials G2_responses2nogo_trials 
% % all subjects data
% G_RT_median=median(G_RT);
% % Subjects_G_RT(s-2,1)=G_RT_median;
% 
% todelete=find(G_RT(:,1)==-1);
% G_RT(todelete, :)=[];
% 
% % plot
% figure; plot(G_RT(:,1), G_RT(:,2),'d');
% filename=strcat('AB', num2str(s), ' GNG');
% title(filename,'FontSize',16,'FontWeight','bold');
% xlabel('Reaction time (s)','FontSize',16,'FontWeight','bold');
% ylabel('Warning to imperative time (s)','FontSize',16,'FontWeight','bold')
% 
% [R, P]=corrcoef(G_RT(:,1), G_RT(:,2));
% % [RHO,PVAL] = corr(D_RT(:,1), D_RT(:,2), 'type', 'Spearman')
% 
% Correlation_RT_Timing_GNG(b)=R(1,2);

end


% cd('M:\ProjectAgingNeuromodulation\AuditoryResearch\BehavioralAnalysis')
% save Correlation_RT_Timing_GNG Correlation_RT_Timing_GNG

% [h p]=ttest(Correlation_RT_Timing_GNG)
