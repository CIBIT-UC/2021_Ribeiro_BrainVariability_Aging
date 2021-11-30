% STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
clear; close all;
younger=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];
x=0; y=0;
task_merged={'W', 'D', 'G'};
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\channel_labels.mat');

% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

for p=8%[younger older]
% clearvars -except t p x y younger older task_merged STUDY CURRENTSTUDY ALLEEG EEG CURRENTSET ...
% CNV_Detection_Young  CNV_Detection_Older  CNV_GNG_Young  CNV_GNG_Older channel_labels

participant=strcat('AB', num2str(p));
directory1=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\EEG\');
% cd(directory1);
% load file in eeglab - detection task
t=2;
% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

filename=strcat(directory1, participant, '_', task_merged{t}, '_Merged_epochs_baselined.set');
EEG = pop_loadset(filename); % data from G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\eeglab_analysis_3.m
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

% EEG.data(channel, time, trials) - to make sure all channels are in the
% right order - might have been changed due to interpolation of bad
% channels - I don't think channels were interpolated in these files -
% CHECK!!!
data_tmp=NaN(length(channel_labels), size(EEG.data, 2), size(EEG.data, 3));
for c=1:size(EEG.chanlocs, 1)
    for cc=1:length(channel_labels)
        if strcmp(EEG.chanlocs(c).labels, channel_labels{cc})
            data_tmp(cc, :, :)=EEG.data(c, :, :);
        end
    end
end


% include only correct go trials and exclude trials after errors!
dir_beh_detection = strcat(directory1,'detection_behavior');
clear misses multiple_responses response2cue slow_responses response2nogo error_trials;
% error_trials - all errors, response2cue, multiple responses, misses, slow responses
load([dir_beh_detection filesep 'misses.mat']);
load([dir_beh_detection filesep 'multiple_responses.mat']);
load([dir_beh_detection filesep 'response2cue.mat']);
load([dir_beh_detection filesep 'slow_responses.mat']);

error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses);

load([dir_beh_detection filesep 'correct_trial.mat']);

% avg erp amplitude in the time window 1000-1500 ms after cue onset (before
% target onset) - epoched data from -.2 up to 6s after cue onset - sampling
% rate 500 Hz
% erp_avg_amplitude = squeeze(mean(data_tmp(:, 351:850, setdiff(correct_trial, error_trials+1)), 2));
erp_avg_amplitude = squeeze(mean(data_tmp(:, 601:850, setdiff(correct_trial, error_trials+1)), 2));

%create variables with standard deviation of data from all subjects
if isempty(setdiff(p, younger))
    x=x+1;
    ERP_StDev_Detection_Young(x, :, :)=std(data_tmp(:, :, setdiff(correct_trial, error_trials+1)), [], 3); % participant x channel x time
    ERP_avg_amp_stdev_young(1, x, :) = std(erp_avg_amplitude, [], 2); % task x participant x channel
    ERP_avg_amp_young(1, x, :) = mean(erp_avg_amplitude, 2); % task x participant x channel
elseif isempty(setdiff(p, older))
    y=y+1;
    ERP_StDev_Detection_Older(y, :, :)=std(data_tmp(:, :, setdiff(correct_trial, error_trials+1)), [], 3); % for each time point
    ERP_avg_amp_stdev_older(1, y, :) = std(erp_avg_amplitude, [], 2); % for avg ERP amplitude
    ERP_avg_amp_older(1, y, :) = mean(erp_avg_amplitude, 2); % task x participant x channel
end

%     % graph of all trials per participant - FCz
%     data_correct = data_tmp(:, :, setdiff(correct_trial, error_trials+1)); x_axis=-0.199:1/500:6;
%     for trl = 1:size(data_correct, 3)
%         clr(trl) = floor(255/size(data_correct, 3)*trl)/255;
%     end
%     color_index = randperm(size(data_correct, 3));
%     figure;
%     for trl = 1:size(data_correct, 3)
%         plot(x_axis, data_correct(16, :, trl), 'color', [clr(color_index(trl)) .5 .5]); hold on
%     end
%     plot(x_axis, mean(data_correct(16, :, :), 3), 'k', 'LineWidth', 2); hold on
%     plot(x_axis, zeros(length(x_axis), 1), '--k')
%     hold off
%     ax = gca; ax.FontSize = 20; ax.FontName = 'Arial'; ax.Color = 'none';
%     axis([-0.2 2 -60 60])
%     xlabel('Time (s)', 'FontSize', 28, 'FontWeight','normal')
%     ylabel('Amplitude /mV', 'FontSize', 28, 'FontWeight','normal')
%     title(['AB', num2str(p), ' - Simple RT'], 'FontSize', 28, 'FontWeight','normal')
%     
%     figure; % plot ERP of participant
%     plot(x_axis, mean(data_correct(16, :, :), 3), 'k', 'LineWidth', 2); hold on
%     plot(x_axis, zeros(length(x_axis), 1), '--k')
%     hold off
%     ax = gca; ax.FontSize = 20; ax.FontName = 'Arial'; ax.Color = 'none';
%     axis([-0.2 2 -20 10])
%     xlabel('Time (s)', 'FontSize', 28, 'FontWeight','normal')
%     ylabel('Amplitude /mV', 'FontSize', 28, 'FontWeight','normal')
%     title(['AB', num2str(p), ' - Simple RT'], 'FontSize', 28, 'FontWeight','normal')
%     
%     %plot 12 random trial to visualize variability
%     figure; trl_numb = randperm(size(data_correct, 3)); count = 0;
%     for trl = 1:12
%         subplot(3, 4, trl);
%         plot(x_axis, data_correct(16, :, trl_numb(trl)), 'color', [0.5430 0 0]); hold on
%         plot(x_axis, zeros(length(x_axis), 1), '--k')
%         hold off; axis([-0.2 2 -40 40])
%     end
%     ax = gca; ax.FontSize = 20; ax.FontName = 'Arial'; ax.Color = 'none';

    
% go/no-go task
%load file in eeglab - go/no-go task
t=3;
% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
% clearvars -except t p x y  participant younger older task_merged directory1 STUDY CURRENTSTUDY ALLEEG EEG CURRENTSET ...
%      CNV_Detection_Young  CNV_Detection_Older  CNV_GNG_Young  CNV_GNG_Older channel_labels
% cd(directory1);
filename=strcat(directory1, participant, '_', task_merged{t}, '_Merged_epochs_baselined.set');
EEG = pop_loadset(filename);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

% EEG.data(channel, time, trials)
data_tmp=NaN(length(channel_labels), size(EEG.data, 2), size(EEG.data, 3));
for c=1:size(EEG.chanlocs, 1)
    for cc=1:length(channel_labels)
        if strcmp(EEG.chanlocs(c).labels, channel_labels{cc})
            data_tmp(cc, :, :)=EEG.data(c, :, :);
        end
    end
end

% include only correct go trials and exclude trials after errors!
clear misses multiple_responses response2cue slow_responses response2nogo error_trials correct_trial;
dir_beh_gng = strcat(directory1,'GNG_behavior');
load([dir_beh_gng filesep 'misses.mat']);
load([dir_beh_gng filesep 'multiple_responses.mat']);
load([dir_beh_gng filesep 'response2cue.mat']);
load([dir_beh_gng filesep 'slow_responses.mat']);
load([dir_beh_gng filesep 'response2nogo.mat']);

% error_trials - all errors, response2cue, multiple responses, misses, slow responses
error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses, response2nogo);

load([dir_beh_gng filesep 'correct_trial.mat']);

% avg erp amplitude in the time window 1000-1500 ms after cue onset (before
% target onset)
erp_avg_amplitude = squeeze(mean(data_tmp(:, 601:850, setdiff(correct_trial, error_trials+1)), 2)); 

    %create variables with mean data from all subjects
    if isempty(setdiff(p, younger))
        ERP_StDev_GNG_Young(x, :, :)=std(data_tmp(:, :, setdiff(correct_trial, error_trials+1)), [], 3);
        ERP_avg_amp_stdev_young(2, x, :) = std(erp_avg_amplitude, [], 2); % task x participant x channel
        ERP_avg_amp_young(2, x, :) = mean(erp_avg_amplitude, 2); % task x participant x channel
    elseif isempty(setdiff(p, older))
        ERP_StDev_GNG_Older(y, :, :)=std(data_tmp(:, :, setdiff(correct_trial, error_trials+1)), [], 3);
        ERP_avg_amp_stdev_older(2, y, :) = std(erp_avg_amplitude, [], 2); % task x participant x channel
        ERP_avg_amp_older(2, y, :) = mean(erp_avg_amplitude, 2); % task x participant x channel
    end
    
    % graph of all trials per participant - FCz
    data_correct = data_tmp(:, :, setdiff(correct_trial, error_trials+1)); x_axis=-0.199:1/500:6;
    for trl = 1:size(data_correct, 3)
        clr(trl) = floor(255/size(data_correct, 3)*trl)/255;
    end
    color_index = randperm(size(data_correct, 3));
    figure;
    for trl = 1:size(data_correct, 3)
%         plot(x_axis, data_correct(16, :, trl), 'color', [clr(color_index(trl)) .5 .5]); hold on
        plot(x_axis, data_correct(16, :, trl), 'color', [.5 clr(color_index(trl)) .5]); hold on
    end
    plot(x_axis, mean(data_correct(16, :, :), 3), 'k', 'LineWidth', 2); hold on 
    plot(x_axis, zeros(length(x_axis), 1), '--k'); hold off
    ax = gca; ax.FontSize = 20; ax.FontName = 'Arial'; ax.Color = 'none';
    axis([-0.2 2 -60 60])
    xlabel('Time (s)', 'FontSize', 28, 'FontWeight','normal')
    ylabel('Amplitude \muV', 'FontSize', 28, 'FontWeight','normal')
    title(['AB', num2str(p), ' - Go-no/go'], 'FontSize', 28, 'FontWeight','normal')
    
    figure; trl_numb = randperm(size(data_correct, 3)); count = 0;
    for trl = 1:12
        subplot(3, 4, trl);
        plot(x_axis, data_correct(16, :, trl_numb(trl)), 'k');  hold on 
        plot(x_axis, zeros(length(x_axis), 1), '--k'); hold off
        axis([0 .5 -60 60])
    end
    trl_numb;

end

%% save variables
save_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\';
% save([save_dir 'ERP_StDev_Detection_Young.mat'], 'ERP_StDev_Detection_Young');
% save([save_dir 'ERP_StDev_Detection_Older.mat'], 'ERP_StDev_Detection_Older');
% save([save_dir 'ERP_StDev_GNG_Young.mat'], 'ERP_StDev_GNG_Young');
% save([save_dir 'ERP_StDev_GNG_Older.mat'], 'ERP_StDev_GNG_Older');

save([save_dir 'ERP_avg_amp_stdev_young.mat'], 'ERP_avg_amp_stdev_young'); % task x participant x channel
save([save_dir 'ERP_avg_amp_stdev_older.mat'], 'ERP_avg_amp_stdev_older');

save([save_dir 'ERP_avg_amp_young.mat'], 'ERP_avg_amp_young');
save([save_dir 'ERP_avg_amp_older.mat'], 'ERP_avg_amp_older');

%% plot data - stdev
save_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\';
load([save_dir 'ERP_StDev_Detection_Young.mat']); load([save_dir 'ERP_StDev_Detection_Older.mat']);
load([save_dir 'ERP_StDev_GNG_Young.mat']); load([save_dir 'ERP_StDev_GNG_Older.mat']);

% FCz = 16
plot_data_mean_std_young_older(squeeze(ERP_StDev_Detection_Young(:, 16,101:3100)), squeeze(ERP_StDev_Detection_Older(:, 16,101:3100)), 'Simple RT', {'Standard deviation' 'ERP'});
plot_data_mean_std_young_older(squeeze(ERP_StDev_GNG_Young(:, 16,101:3100)), squeeze(ERP_StDev_GNG_Older(:, 16,101:3100)), 'Go/no-go', {'Standard deviation' 'ERP'});


%% calculate mean std in epoch and compare across groups - avg variability between 500 and 1500 ms after cue onset
% channel FCz = 16
avg_std_FCz_young(:, 1) = mean(ERP_StDev_Detection_Young(:, 16, 351:850), 3);
avg_std_FCz_older(:, 1) = mean(ERP_StDev_Detection_Older(:, 16, 351:850), 3);
avg_std_FCz_young(:, 2) = mean(ERP_StDev_GNG_Young(:, 16, 351:850), 3);
avg_std_FCz_older(:, 2) = mean(ERP_StDev_GNG_Older(:, 16, 351:850), 3);

% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(avg_std_FCz_young(:, 1), avg_std_FCz_young(:, 2),...
    avg_std_FCz_older(:, 1) , avg_std_FCz_older(:, 2), 'ERP StD (\muV)');

%% standard deviation of avg amplitude of ERP - FCz
save_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\';
load([save_dir 'ERP_avg_amp_stdev_young.mat']);
load([save_dir 'ERP_avg_amp_stdev_older.mat']);
% FCz = 16
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
% plot_all_data_2tasks(squeeze(ERP_avg_amp_stdev_young(1, :, 16)), squeeze(ERP_avg_amp_stdev_young(2, :, 16)),...
%     squeeze(ERP_avg_amp_stdev_older(1, :, 16)), squeeze(ERP_avg_amp_stdev_older(2, :, 16)), 'Amplitude StD (\muV)');
plot_all_data_onetask(squeeze(ERP_avg_amp_stdev_young(2, :, 16)),squeeze(ERP_avg_amp_stdev_older(2, :, 16)), 'Amplitude SD (\muV)')

load([save_dir 'ERP_avg_amp_young.mat']);
load([save_dir 'ERP_avg_amp_older.mat']);
% plot_all_data_2tasks(squeeze(ERP_avg_amp_young(1, :, 16)), squeeze(ERP_avg_amp_young(2, :, 16)),...
%     squeeze(ERP_avg_amp_older(1, :, 16)), squeeze(ERP_avg_amp_older(2, :, 16)), 'Amplitude (\muV)');


plot_all_data_onetask(squeeze(ERP_avg_amp_young(2, :, 16)), squeeze(ERP_avg_amp_older(2, :, 16)), 'Amplitude (\muV)')

%% standard deviation of avg amplitude of ERP
save_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\';
load([save_dir 'ERP_avg_amp_stdev_young.mat']);
load([save_dir 'ERP_avg_amp_stdev_older.mat']);
% FC5 = 13
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(squeeze(ERP_avg_amp_stdev_young(1, :, 13)), squeeze(ERP_avg_amp_stdev_young(2, :, 13)),...
    squeeze(ERP_avg_amp_stdev_older(1, :, 13)), squeeze(ERP_avg_amp_stdev_older(2, :, 13)), 'Amplitude StD (\muV)');

load([save_dir 'ERP_avg_amp_young.mat']);
load([save_dir 'ERP_avg_amp_older.mat']);
plot_all_data_2tasks(squeeze(ERP_avg_amp_young(1, :, 13)), squeeze(ERP_avg_amp_young(2, :, 13)),...
    squeeze(ERP_avg_amp_older(1, :, 13)), squeeze(ERP_avg_amp_older(2, :, 13)), 'Amplitude (\muV)');


%% create excel file for SPSS analyses
T_values(:, 1) = [younger'; older'];
T_values(:, 2) = [ones(size(avg_std_FCz_young, 1), 1); ones(size(avg_std_FCz_older, 1), 1)*2];
T_values(:, 3) = [avg_std_FCz_young(:, 1); avg_std_FCz_older(:, 1)];
T_values(:, 4) = [avg_std_FCz_young(:, 2); avg_std_FCz_older(:, 2)];

T_values(:, 5) = [squeeze(ERP_avg_amp_stdev_young(1, :, 16))'; squeeze(ERP_avg_amp_stdev_older(1, :, 16))'];
T_values(:, 6) = [squeeze(ERP_avg_amp_stdev_young(2, :, 16))'; squeeze(ERP_avg_amp_stdev_older(2, :, 16))'];

T_values(:, 7) = [squeeze(ERP_avg_amp_young(1, :, 16))'; squeeze(ERP_avg_amp_older(1, :, 16))'];
T_values(:, 8) = [squeeze(ERP_avg_amp_young(2, :, 16))'; squeeze(ERP_avg_amp_older(2, :, 16))'];

T = array2table(T_values, ...
    'VariableNames',{'part_id', 'group', 'FCz_erp_std_simpleRT', 'FCz_erp_std_gng', 'FCz_avg_amp_std_simpleRT', 'FCz_avg_amp_std_gng'...
    'FCz_avg_amp_simpleRT', 'FCz_avg_amp_gng'});

filename = 'erp_standard_deviation_500ms.xlsx';
writetable(T,filename,'Sheet',1,'Range','A1')

% merge with excel file with pre-stimulus eeg signal variability and
% spectral parameters for stats in spss calculated in:
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\foof_analysis_power_spectrum_FCz_thresh1.m
% and
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\pre_stimulus_variability_3sec.m

% for statistical control use zero power for subjects where alpha oscillations were not
% detected


%% plot data
save_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\';
load([save_dir 'ERP_StDev_Detection_Young.mat']); load([save_dir 'ERP_StDev_Detection_Older.mat']);
load([save_dir 'ERP_StDev_GNG_Young.mat']); load([save_dir 'ERP_StDev_GNG_Older.mat']);
Pz = 43;
plot_data_mean_std_young_older(squeeze(ERP_StDev_Detection_Young(:, Pz,101:3100)), squeeze(ERP_StDev_Detection_Older(:, Pz,101:3100)), 'Simple RT', {'Standard deviation' 'ERP'});
plot_data_mean_std_young_older(squeeze(ERP_StDev_GNG_Young(:, Pz,101:3100)), squeeze(ERP_StDev_GNG_Older(:, Pz,101:3100)), 'Go/no-go', {'Standard deviation' 'ERP'});



%% create graphs

function plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)

    % plot data for young group - simple RT and go/nogo task
    figure; box on; hold on
    
    % plot data for group 1
    yMean1=nanmean(data_grp1_task1); yMean2=nanmean(data_grp1_task2);
    y_se1 = std(data_grp1_task1)/sqrt(length(data_grp1_task1)); y_se2 = std(data_grp1_task2)/sqrt(length(data_grp1_task2));
    
    %plot the mean+-SEM box
    %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
    %   Specify pos as a four-element vector of the form [x y w h] in data
    %   units. The x and y elements determine the location and the w and h
    %   elements determine the size. The function plots into the current axes
    %   without clearing existing content from the axes.
    rectangle('Position',[1-0.3,yMean1-y_se1, 0.6, 2*y_se1 ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    rectangle('Position',[2-0.3,yMean2-y_se2, 0.6, 2*y_se2 ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    
    for y=1:length(data_grp1_task1)
        plot([1 2]+rand*0.2-0.1, [data_grp1_task1(y) data_grp1_task2(y)] ,'-o', 'color', [.8 .8 .8], ...
            'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end

    %plot the mean line
    plot([1 2 ],[yMean1 yMean2] ,'Color','k','LineWidth',1.5);
    plot([1-0.3 1+0.3],[yMean1 yMean1] ,'Color','k','LineWidth',5);
    plot([2-0.3 2+0.3],[yMean2 yMean2] ,'Color','k','LineWidth',5);
    
    
    %% group 2
        yMean1=nanmean(data_grp2_task1); yMean2=nanmean(data_grp2_task2);
    y_se1 = std(data_grp2_task1)/sqrt(length(data_grp2_task1)); y_se2 = std(data_grp2_task2)/sqrt(length(data_grp2_task2));
    
    %plot the mean+-SEM box
    %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
    %   Specify pos as a four-element vector of the form [x y w h] in data
    %   units. The x and y elements determine the location and the w and h
    %   elements determine the size. The function plots into the current axes
    %   without clearing existing content from the axes.
    rectangle('Position',[4-0.3,yMean1-y_se1, 0.6, 2*y_se1 ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    rectangle('Position',[5-0.3,yMean2-y_se2, 0.6, 2*y_se2 ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    
    for y=1:length(data_grp2_task1)
        plot([4 5]+rand*0.2-0.1, [data_grp2_task1(y) data_grp2_task2(y)] ,'-o', 'color', [1 .8 .8], ...
            'MarkerFaceColor',[1 .8 .8], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;  
    end

     %plot the mean line
    plot([4 5],[yMean1 yMean2] ,'Color','k','LineWidth',1.5);
    plot([4-0.3 4+0.3],[yMean1 yMean1] ,'Color','k','LineWidth',5);
    plot([5-0.3 5+0.3],[yMean2 yMean2] ,'Color','k','LineWidth',5);
    

    % axes('XColor','none');
    hold off;
    axis([0 6 -inf inf]);
    ax = gca;
    c = ax.Color;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel=[1 2 1 2];
    xticks([1 2 4 5])
    ylabel(y_label_text, 'FontSize', 32, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end


function plot_data_mean_std_young_older(data_young, data_older, graph_title, y_axis_label)

mean_data_young=mean(data_young, 1);
se_data_young=std(data_young, 0, 1)/sqrt(size(data_young, 1));

mean_data_older=mean(data_older, 1);
se_data_older=std(data_older, 0, 1)/sqrt(size(data_older, 1));


% color gng young = [0 .5 0]
% color simple RT young =   [1 .5 0]
% color gng older = [0 0 .5]
% color simple RT older =  [.75 0.25 0]

x_axis=0.001:1/500:6;
figure;
% % plot a line at zero
% plot([0 0],[0 11], '--', 'color', [0 0 0]);
% hold on
% detection data
plot(x_axis, mean_data_young,'color', 'k', 'LineWidth',1.5)
hold on
jbfill(x_axis,mean_data_young+se_data_young, mean_data_young-se_data_young, 'k','k', 0.1)
% GNG data
hold on
plot(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
hold on
jbfill(x_axis, mean_data_older+se_data_older, mean_data_older-se_data_older,'r','r', 0.1)
hold off
ax = gca;
% c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 20;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-inf inf -inf inf])
xlabel('Time (s)', 'FontSize', 28, 'FontWeight','normal')
ylabel(y_axis_label, 'FontSize', 28, 'FontWeight','normal')
title(graph_title, 'FontSize', 28, 'FontWeight','normal')

end

% plot data from one task only
function plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)

    % plot data for young group - go/nogo task
    figure; box on; hold on
    
    % plot data for group 1
    yMean1=nanmean(data_grp1_task1);
    y_se1 = std(data_grp1_task1)/sqrt(length(data_grp1_task1));
    
    %plot the mean+-SEM box
    %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
    %   Specify pos as a four-element vector of the form [x y w h] in data
    %   units. The x and y elements determine the location and the w and h
    %   elements determine the size. The function plots into the current axes
    %   without clearing existing content from the axes.
    rectangle('Position',[1-0.3,yMean1-y_se1, 0.6, 2*y_se1 ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    
    for y=1:length(data_grp1_task1)
        plot(1+rand*0.2-0.1, data_grp1_task1(y),'-o', 'color', [.8 .8 .8], ...
            'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5);
        hold on;
    end

    %plot the mean line
    plot(1,yMean1,'Color','k','LineWidth',1.5);
    plot([1-0.3 1+0.3],[yMean1 yMean1] ,'Color','k','LineWidth',5);
    
    
    %% group 2
        yMean1=nanmean(data_grp2_task1);
    y_se1 = std(data_grp2_task1)/sqrt(length(data_grp2_task1));
    
    %plot the mean+-SEM box
    %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
    %   Specify pos as a four-element vector of the form [x y w h] in data
    %   units. The x and y elements determine the location and the w and h
    %   elements determine the size. The function plots into the current axes
    %   without clearing existing content from the axes.
    rectangle('Position',[2-0.3,yMean1-y_se1, 0.6, 2*y_se1 ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    
    for y=1:length(data_grp2_task1)
        plot(2+rand*0.2-0.1, data_grp2_task1(y),'-o', 'color', [1 .5 .5], ...
            'MarkerFaceColor',[1 .5 .5], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5);
        hold on;  
    end

     %plot the mean line
    plot(2,yMean1 ,'Color','k','LineWidth',1.5);
    plot([2-0.3 2+0.3],[yMean1 yMean1] ,'Color','k','LineWidth',5);

    % axes('XColor','none');
    hold off;
    axis([0 3 -inf inf]);
    ax = gca;
    c = ax.Color;
    ax.YAxis.FontSize = 24;
    ax.XAxis.FontSize = 32;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel= {'Young' 'Older'};
    xticks([1 2])
    ylabel(y_label_text, 'FontSize', 32, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end