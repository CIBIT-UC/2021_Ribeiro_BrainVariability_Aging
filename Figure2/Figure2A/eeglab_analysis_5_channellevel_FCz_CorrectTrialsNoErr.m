% STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
clear; close all;
younger=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];
x=0; y=0; z=0; w=0;
task_merged={'W', 'D', 'G'};

for p = [younger older]

participant=strcat('AB', num2str(p));
directory1=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\EEG\');
directory=strcat('G:\\ProjectAgingNeuromodulation\\AuditoryResearch\\AuditoryTask_EyeTracking_EEG_lab94\\', participant, '\\EEG\\');

cd(directory1);

% creating variable with FCz data
% Maria Ribeiro - May2016
clearvars -except t p x y z w younger older participant task_merged directory1 STUDY CURRENTSTUDY ALLEEG EEG CURRENTSET ...
    FCz_Detection_CorrectTrialsNoErr_YoungGroup FCz_Detection_CorrectTrialsNoErr_OlderGroup FCz_GNG_CorrectTrialsNoErr_YoungGroup FCz_GNG_CorrectTrialsNoErr_OlderGroup

% % clear eeglab
% STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

%load file in eeglab - detection task
t=2;
clearvars -except t p x y z w  participant younger older task_merged directory1 STUDY CURRENTSTUDY ALLEEG EEG CURRENTSET ...
    FCz_Detection_CorrectTrialsNoErr_YoungGroup FCz_Detection_CorrectTrialsNoErr_OlderGroup FCz_GNG_CorrectTrialsNoErr_YoungGroup FCz_GNG_CorrectTrialsNoErr_OlderGroup
% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

filename=strcat(participant, '_', task_merged{t}, '_Merged_epochs_baselined.set');
EEG = pop_loadset(filename);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

% FCz data EEG.data(channel, time, trials)
% find FCz 
for c=1:size(EEG.chanlocs, 1)
    if strcmp(EEG.chanlocs(c).labels, 'FCz')
        a=c;
    end
end

% FCz data EEG.data(channel, time, trials)
FCz_Detection=squeeze(EEG.data(a,:,:));

% include only correct go trials and exclude trials after errors
% load behavioural parameters calculated in eeglab_analysis_4_behavior.m
cd(strcat(pwd,'/detection_behavior'));
clear misses multiple_responses response2cue slow_responses response2nogo error_trials;
% error_trials - all errors, response2cue, multiple responses, misses, slow responses

load misses.mat;
load multiple_responses.mat;
load response2cue.mat;
load slow_responses.mat;

error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses);

load correct_trial.mat

% channel of interest FCz - only correct go trials, excluding trials
% after error

FCz_Detection_CorrectTrialsNoErr=FCz_Detection( :, setdiff(correct_trial, error_trials+1));
cd(directory1);
save FCz_Detection_CorrectTrialsNoErr FCz_Detection_CorrectTrialsNoErr

load FCz_Detection_CorrectTrialsNoErr

%create variables with mean data from all subjects
if isempty(setdiff(p, younger))
    x=x+1;
    FCz_Detection_CorrectTrialsNoErr_YoungGroup(x, :)=mean(FCz_Detection_CorrectTrialsNoErr, 2);
elseif isempty(setdiff(p, older))
    y=y+1;
    FCz_Detection_CorrectTrialsNoErr_OlderGroup(y, :)=mean(FCz_Detection_CorrectTrialsNoErr, 2);
end


% go/no-go task
%load file in eeglab - go/no-go task
t=3;
% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
clearvars -except t p a x y z w  participant younger older task_merged directory1 STUDY CURRENTSTUDY ALLEEG EEG CURRENTSET ...
    FCz_Detection_CorrectTrialsNoErr_YoungGroup FCz_Detection_CorrectTrialsNoErr_OlderGroup FCz_GNG_CorrectTrialsNoErr_YoungGroup FCz_GNG_CorrectTrialsNoErr_OlderGroup

filename=strcat(participant, '_', task_merged{t}, '_Merged_epochs_baselined.set');
EEG = pop_loadset(filename);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

% FCz data EEG.data(channel, time, trials)
FCz_GNG=squeeze(EEG.data(a,:,:));

% include only correct go trials and exclude trials after errors
% load behavioural parameters calculated in eeglab_analysis_4_behavior.m
clear misses multiple_responses response2cue slow_responses response2nogo error_trials correct_trial;
cd(strcat(pwd,'/GNG_behavior'));
load misses.mat;
load multiple_responses.mat;
load response2cue.mat;
load slow_responses.mat;
load response2nogo.mat

% error_trials - all errors, response2cue, multiple responses, misses, slow responses
error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses, response2nogo);

load correct_trial.mat

% channel of interest FCz - only correct go trials, excluding trials
% after error

FCz_GNG_CorrectTrialsNoErr=FCz_GNG( :, setdiff(correct_trial, error_trials+1));
cd(directory1);
save FCz_GNG_CorrectTrialsNoErr FCz_GNG_CorrectTrialsNoErr

load FCz_GNG_CorrectTrialsNoErr

%create variables with mean data from all subjects
if isempty(setdiff(p, younger))
    z=z+1;
    FCz_GNG_CorrectTrialsNoErr_YoungGroup(z, :)=mean(FCz_GNG_CorrectTrialsNoErr, 2);
elseif isempty(setdiff(p, older))
    w=w+1;
    FCz_GNG_CorrectTrialsNoErr_OlderGroup(w, :)=mean(FCz_GNG_CorrectTrialsNoErr, 2);
end

end

%% save group variables

cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\FCz_group_data');

save('FCz_Detection_CorrectTrialsNoErr_OlderGroup', 'FCz_Detection_CorrectTrialsNoErr_OlderGroup')
save('FCz_Detection_CorrectTrialsNoErr_YoungGroup', 'FCz_Detection_CorrectTrialsNoErr_YoungGroup')
save('FCz_GNG_CorrectTrialsNoErr_OlderGroup', 'FCz_GNG_CorrectTrialsNoErr_OlderGroup')
save('FCz_GNG_CorrectTrialsNoErr_YoungGroup', 'FCz_GNG_CorrectTrialsNoErr_YoungGroup')

%% plot graphs

% figure with standard errors - detection older vs younger
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\FCz_group_data');
load FCz_Detection_CorrectTrialsNoErr_YoungGroup;
load FCz_Detection_CorrectTrialsNoErr_OlderGroup;

MeanYoung=mean(FCz_Detection_CorrectTrialsNoErr_YoungGroup, 1);
SEYoung=std(FCz_Detection_CorrectTrialsNoErr_YoungGroup,0, 1)/sqrt(size(FCz_Detection_CorrectTrialsNoErr_OlderGroup, 1));
MeanOlder=mean(FCz_Detection_CorrectTrialsNoErr_OlderGroup, 1);
SEOlder=std(FCz_Detection_CorrectTrialsNoErr_OlderGroup,0, 1)/sqrt(size(FCz_Detection_CorrectTrialsNoErr_OlderGroup, 1));

x_axis=[];
x_axis=-200:1000/500:5999;
figure;
for x = 1001:2:1500
    plot([x x],[-13.8 1.8], 'color', [.8 .8 .8] ); hold on
end
hold on
% detection data
plot(x_axis, MeanYoung,'color', 'k', 'LineWidth',1.5)
hold on
jbfill(x_axis,MeanYoung+SEYoung,MeanYoung-SEYoung, 'k','k', 0.1)
% GNG data
hold on
plot(x_axis, MeanOlder,'color', 'r','LineStyle', '--','LineWidth',1.5)
hold on
jbfill(x_axis,MeanOlder+SEOlder,MeanOlder-SEOlder, 'r','r', 0.1)
hold on
plot(x_axis, zeros(length(x_axis), 1),'color', [0.5 0.5 0.5],'LineStyle', '--', 'LineWidth',1.5)
hold off
box off
ax = gca;
ax.LineWidth = 2.5; 
% legend('Detection', 'GNG')
ax.FontSize = 28;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-200 3000 -14 2])
xticks([0 1000 2000])
xticklabels({'0','1000','2000'})
xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
title('FCz - simple RT' , 'FontSize', 32, 'FontWeight','normal')


%% figure with standard errors - gng older vs younger

cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\FCz_group_data');
load FCz_GNG_CorrectTrialsNoErr_YoungGroup;
load FCz_GNG_CorrectTrialsNoErr_OlderGroup;

MeanYoung=mean(FCz_GNG_CorrectTrialsNoErr_YoungGroup, 1);
SEYoung=std(FCz_GNG_CorrectTrialsNoErr_YoungGroup,0, 1)/sqrt(size(FCz_GNG_CorrectTrialsNoErr_OlderGroup, 1));
MeanOlder=mean(FCz_GNG_CorrectTrialsNoErr_OlderGroup, 1);
SEOlder=std(FCz_GNG_CorrectTrialsNoErr_OlderGroup,0, 1)/sqrt(size(FCz_GNG_CorrectTrialsNoErr_OlderGroup, 1));

x_axis=-200:1000/500:5999;
x_axis = x_axis/1000; % in seconds
figure;
for x = 1.001:.002:1.500
    plot([x x],[-13.8 1.8], 'color', [.8 .8 .8] ); hold on
end
hold on
% detection data
plot(x_axis, MeanYoung,'color', 'k',  'LineWidth',1.5)
hold on
jbfill(x_axis,MeanYoung+SEYoung,MeanYoung-SEYoung, 'k','k', 0.1)
% GNG data
hold on
plot(x_axis, MeanOlder,'color', 'r','LineStyle', '--','LineWidth',1.5)
hold on
jbfill(x_axis,MeanOlder+SEOlder,MeanOlder-SEOlder,'r','r', 0.1)
hold on
plot(x_axis, zeros(length(x_axis), 1),'color', [0.5 0.5 0.5],'LineStyle', '--', 'LineWidth',1.5)
hold off
box off
ax = gca;
ax.LineWidth = 2.5; 
% legend('GNG', 'GNG')
ax.FontSize = 28;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-.200 3 -14 2])
xticks([0 1 2]);
xticklabels({'0','1','2'})
% axis([0 500 -14 2])
xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
title('FCz' , 'FontSize', 32, 'FontWeight','normal')

%% detection vs gng younger
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\FCz_group_data')
load FCz_GNG_CorrectTrialsNoErr_YoungGroup;
load FCz_GNG_CorrectTrialsNoErr_OlderGroup;
load FCz_Detection_CorrectTrialsNoErr_YoungGroup;
load FCz_Detection_CorrectTrialsNoErr_OlderGroup;

% colors ok for colorblind
% detection = [.9 .6 0]
% gng = [.35 .7 .9]

MeanYoung_detection=mean(FCz_Detection_CorrectTrialsNoErr_YoungGroup, 1);
SEYoung_detection=std(FCz_Detection_CorrectTrialsNoErr_YoungGroup,0, 1)/sqrt(size(FCz_Detection_CorrectTrialsNoErr_YoungGroup, 1));

MeanYoung_gng=mean(FCz_GNG_CorrectTrialsNoErr_YoungGroup, 1);
SEYoung_gng=std(FCz_GNG_CorrectTrialsNoErr_YoungGroup,0, 1)/sqrt(size(FCz_GNG_CorrectTrialsNoErr_YoungGroup, 1));


x_axis=[];
x_axis=-200:1000/500:5999;
figure;
% plot line at zero
plot(x_axis, zeros(length(x_axis), 1),'color', [0.5 0.5 0.5],'LineStyle', '--', 'LineWidth',1.5)
hold on
% detection data
plot(x_axis, MeanYoung_detection,'color', [.9 .6 0], 'LineWidth',1.5)
hold on
jbfill(x_axis,MeanYoung_detection+SEYoung_detection,MeanYoung_detection-SEYoung_detection,  [.9 .6 0], [.9 .6 0], 0.1)
% GNG data
hold on
plot(x_axis, MeanYoung_gng,'color', [.35 .7 .9],'LineStyle', '--', 'LineWidth',1.5)
hold on
jbfill(x_axis,MeanYoung_gng+SEYoung_gng,MeanYoung_gng-SEYoung_gng,[.35 .7 .9],[.35 .7 .9], 0.1)
hold off
ax = gca;
c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 28;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-200 3000 -14 2])
xticks([0 1000 2000])
xticklabels({'0','1000','2000'})
xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','bold')
ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','bold')
title('FCz - young group' , 'FontSize', 32, 'FontWeight','bold')

%% older detection vs gng
MeanOlder_detection=mean(FCz_Detection_CorrectTrialsNoErr_OlderGroup, 1);
SEOlder_detection=std(FCz_Detection_CorrectTrialsNoErr_OlderGroup,0, 1)/sqrt(size(FCz_Detection_CorrectTrialsNoErr_OlderGroup, 1));

MeanOlder_gng=mean(FCz_GNG_CorrectTrialsNoErr_OlderGroup, 1);
SEOlder_gng=std(FCz_GNG_CorrectTrialsNoErr_OlderGroup,0, 1)/sqrt(size(FCz_GNG_CorrectTrialsNoErr_OlderGroup, 1));


x_axis=[];
x_axis=-200:1000/500:5999;
figure;
% plot line at zero
plot(x_axis, zeros(length(x_axis), 1),'color', [0.5 0.5 0.5],'LineStyle', '--', 'LineWidth',1.5)
hold on
% detection data
plot(x_axis, MeanOlder_detection,'color',[.9 .6 0], 'LineWidth',1.5)
hold on
jbfill(x_axis,MeanOlder_detection+SEOlder_detection,MeanOlder_detection-SEOlder_detection, [.9 .6 0],[.9 .6 0], 0.1)
% GNG data
hold on
plot(x_axis, MeanOlder_gng,'color',[.35 .7 .9],'LineStyle', '--', 'LineWidth',1.5)
hold on
jbfill(x_axis,MeanOlder_gng+SEOlder_gng,MeanOlder_gng-SEOlder_gng,[.35 .7 .9],[.35 .7 .9], 0.1)
hold off
ax = gca;
c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 20;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-200 3000 -14 2])
xticks([0 1000 2000])
xticklabels({'0','1000','2000'})
xlabel('Time (ms)', 'FontSize', 28, 'FontWeight','bold')
ylabel('Amplitude (\muV)', 'FontSize', 28, 'FontWeight','bold')
title('FCz - older group' , 'FontSize', 28, 'FontWeight','bold')

%% individual ERPs
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\FCz_group_data')
load FCz_GNG_CorrectTrialsNoErr_YoungGroup;
load FCz_GNG_CorrectTrialsNoErr_OlderGroup;
load FCz_Detection_CorrectTrialsNoErr_YoungGroup;
load FCz_Detection_CorrectTrialsNoErr_OlderGroup;
% older detection vs gng
MeanOlder_detection=mean(FCz_Detection_CorrectTrialsNoErr_OlderGroup, 1);

x_axis=[];
x_axis=-200:1000/500:5999;
figure;

% detection data
for i=1:size(FCz_Detection_CorrectTrialsNoErr_OlderGroup, 1)
    plot(x_axis, FCz_Detection_CorrectTrialsNoErr_OlderGroup(i, :), 'color', [0.75 0.75 0.75],'LineWidth',.5);
    hold on
end
plot(x_axis, MeanOlder_detection, 'color', 'k','LineWidth',2);
hold on
% plot a line at zero
plot([0 0],[-30 12], '--', 'color','k','LineWidth',.5);
hold off
ax = gca;
c = ax.Color;
ax.FontSize = 20;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-200 3000 -inf inf])
% xticks([0 500 1000])
% xticklabels({'0', '500', '1000'})
xlabel('Time (ms)', 'FontSize', 28, 'FontWeight','bold')
ylabel('Amplitude (\muV)', 'FontSize', 28, 'FontWeight','bold')
title('FCz - older group - simple RT' , 'FontSize', 24, 'FontWeight','bold')

%% individual ERPs
% older gng
MeanOlder_gng=mean(FCz_GNG_CorrectTrialsNoErr_OlderGroup, 1);

x_axis=[];
x_axis=-200:1000/500:5999;
figure;
% detection data
for i=1:size(FCz_GNG_CorrectTrialsNoErr_OlderGroup, 1)
    plot(x_axis, FCz_GNG_CorrectTrialsNoErr_OlderGroup(i, :), 'color', [0.75 0.75 0.75],'LineWidth',.5);
    hold on
end
plot(x_axis, MeanOlder_gng, 'color', 'k','LineWidth',2);
hold on
% plot a line at zero
plot([0 0],[-25 10], '--', 'color','k','LineWidth',.5);
hold off
ax = gca;
c = ax.Color;
ax.FontSize = 20;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-200 3000 -inf inf])
xlabel('Time (ms)', 'FontSize', 28, 'FontWeight','bold')
ylabel('Amplitude (\muV)', 'FontSize', 28, 'FontWeight','bold')
title('FCz - older group - go/no-go' , 'FontSize', 24, 'FontWeight','bold')

%% individual ERPs
% young detection vs gng
MeanYoung_detection=mean(FCz_Detection_CorrectTrialsNoErr_YoungGroup, 1);

x_axis=[];
x_axis=-200:1000/500:5999;
figure;
% detection data
for i=1:size(FCz_Detection_CorrectTrialsNoErr_YoungGroup, 1)
    plot(x_axis, FCz_Detection_CorrectTrialsNoErr_YoungGroup(i, :), 'color', [0.75 0.75 0.75],'LineWidth',.5);
    hold on
end
plot(x_axis, MeanYoung_detection, 'color', 'k','LineWidth',2);
hold on
% plot a line at zero
plot([0 0],[-28 10], '--', 'color','k','LineWidth',.5);
hold off
ax = gca;
c = ax.Color;
ax.FontSize = 20;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-200 3000 -inf inf])
xlabel('Time (ms)', 'FontSize', 28, 'FontWeight','bold')
ylabel('Amplitude (\muV)', 'FontSize', 28, 'FontWeight','bold')
title('FCz - young group - simple RT' , 'FontSize', 24, 'FontWeight','bold')

%% individual ERPs
% young gng
MeanYoung_gng=mean(FCz_GNG_CorrectTrialsNoErr_YoungGroup, 1);

x_axis=[];
x_axis=-200:1000/500:5999;
figure;
% detection data 
for i=1:size(FCz_GNG_CorrectTrialsNoErr_YoungGroup, 1)
    plot(x_axis, FCz_GNG_CorrectTrialsNoErr_YoungGroup(i, :), 'color', [0.75 0.75 0.75],'LineWidth',.5);
    hold on
end
plot(x_axis, MeanYoung_gng, 'color', 'k','LineWidth',2);
hold on
% plot a line at zero
plot([0 0],[-28 15], '--', 'color','k','LineWidth',.5);
hold off
ax = gca;
c = ax.Color;
ax.FontSize = 20;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-200 3000 -inf inf])
xlabel('Time (ms)', 'FontSize', 28, 'FontWeight','bold')
ylabel('Amplitude (\muV)', 'FontSize', 28, 'FontWeight','bold')
title('FCz - young group - go/no-go' , 'FontSize', 24, 'FontWeight','bold')