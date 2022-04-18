% Calculate group RT and accuracy
% Maria Ribeiro - June 2016

% correct reaction time
% median, mean and standard deviation

clear
young=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];%33 behavior data not saved estimated from EEG!
older=[7	8	11	12	14	17	19	20	21	22	23	32	35	37	38  41	43	47	48	49	52	55	57	58	60	61	63	64	65	67	69	70	71	73	75	77	79	83  86]; % 19 and 38 behavior data not saved estimated from EEG

% younger group analysis
b=0;
for s=young
b=b+1;
clearvars -except young older b YoungerGroup_G_RT YoungerGroup_D_RT YoungerGroup_D_errors YoungerGroup_G_errors task s

cd(strcat('E:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(s),'\BehavioralAnalysis'));
%% REACTION TIME ANALYSIS
% load reaction times
load D_RT_WithSlow.mat
YoungerGroup_D_RT(b, 1)=median(D_RT(:,1));
YoungerGroup_D_RT(b, 2)=mean(D_RT(:,1));
YoungerGroup_D_RT(b, 3)=std(D_RT(:,1));

load G_RT_WithSlow.mat
YoungerGroup_G_RT(b, 1)=median(G_RT(:,1));
YoungerGroup_G_RT(b, 2)=mean(G_RT(:,1));
YoungerGroup_G_RT(b, 3)=std(G_RT(:,1));

%% ACCURACY ANALYSIS
% detection task
% YoungerGroup_D_errors == missed responses, responses2cue, slow responses, responses2blank

load D1_MissedResponses.mat; load D2_MissedResponses.mat;
YoungerGroup_D_errors(b, 1)=length(D1_MissedResponses)+length(D2_MissedResponses);

load D1_responses2cue_trials.mat; load D2_responses2cue_trials.mat;
YoungerGroup_D_errors(b, 2)=length(D1_responses2cue_trials)+length(D2_responses2cue_trials);

load D1_slow_responses.mat; load D2_slow_responses.mat;
YoungerGroup_D_errors(b, 3)=length(D1_slow_responses)+length(D2_slow_responses);

load D1_responses2blanks_trials.mat; load D2_responses2blanks_trials.mat;
YoungerGroup_D_errors(b, 4)=length(D1_responses2blanks_trials)+length(D2_responses2blanks_trials);

% GNG TASK
% YoungerGroup_G_errors == missed responses-2, responses2cue, slow responses-2, responses2blanks-2, no-go
% errors
load G1_MissedResponses.mat; load G2_MissedResponses.mat;
YoungerGroup_G_errors(b, 1)=length(G1_MissedResponses)+length(G2_MissedResponses);

load G1_responses2cue_trials.mat; load G2_responses2cue_trials.mat;
YoungerGroup_G_errors(b, 2)=length(G1_responses2cue_trials)+length(G2_responses2cue_trials);

load G1_slow_responses.mat; load G2_slow_responses.mat;
YoungerGroup_G_errors(b, 3)=length(G1_slow_responses)+length(G2_slow_responses);

load G1_responses2blanks_trials.mat; load G2_responses2blanks_trials.mat;
YoungerGroup_G_errors(b, 4)=length(G1_responses2blanks_trials)+length(G2_responses2blanks_trials);

load G1_responses2nogo_trials.mat; load G2_responses2nogo_trials.mat;
YoungerGroup_G_errors(b, 5)=length(G1_responses2nogo_trials)+length(G2_responses2nogo_trials);

end

% saving results
cd E:\ProjectAgingNeuromodulation\AuditoryResearch\BehavioralAnalysis

save YoungerGroup_D_RT YoungerGroup_D_RT
save YoungerGroup_G_RT YoungerGroup_G_RT
save YoungerGroup_D_errors YoungerGroup_D_errors
save YoungerGroup_G_errors YoungerGroup_G_errors

%% Older group analysis
b=0;
for s=older;
b=b+1;
clearvars -except young older b OlderGroup_G_RT OlderGroup_D_RT OlderGroup_D_errors OlderGroup_G_errors task s

cd(strcat('E:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(s),'\BehavioralAnalysis'));
%% REACTION TIME ANALYSIS
% load reaction times
load D_RT_WithSlow.mat
OlderGroup_D_RT(b, 1)=median(D_RT(:,1));
OlderGroup_D_RT(b, 2)=mean(D_RT(:,1));
OlderGroup_D_RT(b, 3)=std(D_RT(:,1));

load G_RT_WithSlow.mat
OlderGroup_G_RT(b, 1)=median(G_RT(:,1));
OlderGroup_G_RT(b, 2)=mean(G_RT(:,1));
OlderGroup_G_RT(b, 3)=std(G_RT(:,1));

%% ACCURACY ANALYSIS
% detection task
% OlderGroup_D_errors == missed responses, responses2cue, slow responses, responses2blank

load D1_MissedResponses.mat; load D2_MissedResponses.mat;
OlderGroup_D_errors(b, 1)=length(D1_MissedResponses)+length(D2_MissedResponses);

load D1_responses2cue_trials.mat; load D2_responses2cue_trials.mat;
OlderGroup_D_errors(b, 2)=length(D1_responses2cue_trials)+length(D2_responses2cue_trials);

load D1_slow_responses.mat; load D2_slow_responses.mat;
OlderGroup_D_errors(b, 3)=length(D1_slow_responses)+length(D2_slow_responses);

load D1_responses2blanks_trials.mat; load D2_responses2blanks_trials.mat;
OlderGroup_D_errors(b, 4)=length(D1_responses2blanks_trials)+length(D2_responses2blanks_trials);

% GNG TASK
% OlderGroup_G_errors == missed responses-2, responses2cue, slow responses-2, responses2blanks-2, no-go
% errors
load G1_MissedResponses.mat; load G2_MissedResponses.mat;
OlderGroup_G_errors(b, 1)=length(G1_MissedResponses)+length(G2_MissedResponses);

load G1_responses2cue_trials.mat; load G2_responses2cue_trials.mat;
OlderGroup_G_errors(b, 2)=length(G1_responses2cue_trials)+length(G2_responses2cue_trials);

load G1_slow_responses.mat; load G2_slow_responses.mat;
OlderGroup_G_errors(b, 3)=length(G1_slow_responses)+length(G2_slow_responses);

load G1_responses2blanks_trials.mat; load G2_responses2blanks_trials.mat;
OlderGroup_G_errors(b, 4)=length(G1_responses2blanks_trials)+length(G2_responses2blanks_trials);

load G1_responses2nogo_trials.mat; load G2_responses2nogo_trials.mat;
OlderGroup_G_errors(b, 5)=length(G1_responses2nogo_trials)+length(G2_responses2nogo_trials);

end

%% saving results
cd E:\ProjectAgingNeuromodulation\AuditoryResearch\BehavioralAnalysis

save OlderGroup_D_RT OlderGroup_D_RT
save OlderGroup_G_RT OlderGroup_G_RT
save OlderGroup_D_errors OlderGroup_D_errors
save OlderGroup_G_errors OlderGroup_G_errors


% %% data with participant ID as first column
% % make excel file with data
% excel_text={'ParticipantID'	'Age'	'Sex'	'Group'	'Task Order'	'DetectionMedianRT'	'GNGMedianRT'	'DifRT'	'DetectionMeanRT'	'GNGMeanRT'	'DetectionStdRT' 'GNGStdRT'	'D misses'	'D responses2cue	D slow responses	D responses2blank	G misses	G responses2cue	G slow responses	G responses2blank	no-go errors
% 
% Part_Older_D_RT(:,1)=older';
% Part_Older_D_RT(:,2:4)=OlderGroup_D_RT;
% 
% Part_Older_G_RT(:,1)=older';
% Part_Older_G_RT(:,2:4)=OlderGroup_G_RT;
% 
% load YoungerGroup_D_RT 
% load YoungerGroup_G_RT
% 
% Part_Younger_D_RT(:,1)=young';
% Part_Younger_D_RT(:,2:4)=YoungerGroup_D_RT;
% 
% Part_Younger_G_RT(:,1)=young';
% Part_Younger_G_RT(:,2:4)=YoungerGroup_G_RT;

