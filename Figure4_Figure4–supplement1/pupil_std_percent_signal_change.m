% across trials standard deviation of amplitude of cue-locked pupil
% response during the preparatory period before target onset
clear; close all;
younger=[4 9 10 13 15 16 25 26 28 31 33 34 36 42 44 45 46 50 51 53 54 56 59 62 66 68 72 74 76 78 80 81 82 84 85];
older=[7 8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77  79 83   86];
task={'D1', 'D2', 'G1', 'G2'};
group = {younger older};
std_pupil = cell(2, 5); std_pupil_avg_response = cell(2, 2); pupil_avg_response = cell(2, 2);
pupil_cuelocked_response = cell(2, 2);
% pupil_percent_signal_change = cell(2, 5); 
% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for grp = 1:2
    count = 0;
    for p = group{grp}
        count = count + 1;
        pupil_directory=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\');
        load([pupil_directory, 'RejectEpochsCorrect']); % data from G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\pupil_analysis_3b_RejectedEpochs_SaveRunData.m
        
%         if intersect(p, [8, 11, 14, 15]) % in these participants use second repetition of recordign due to technical problems
%             task = {'W2', 'D1', 'D2', 'G1', 'G2'};
%         else
%             task={'W1', 'D1', 'D2', 'G1', 'G2'};
%         end
        
        clear pupil_avg_response_tmp pupil_response
        RT = cell(4, 1);
        for t = 1:4 % task runs excluding first run = passive task
        
            % calculate pupil go trials data
            filename=strcat(pupil_directory, 'AB', num2str(p), '_', task(t), '_PupilData.mat');
            load(filename{1})
            % use smoothed pupil diameter data with blink artefacts interpolated
            % measured in horizontal direction (x)
            pupil_data_percent_signal_change = PupilData(:, 3)/mean(PupilData(:, 3)) - 1;
            
%           figure; plot(pupil_data_percent_signal_change);
%             pupil_percent_signal_change{grp, t} = [pupil_percent_signal_change{grp, t}; pupil_data_percent_signal_change(1:50000)'];
            std_pupil{grp, t} = [std_pupil{grp, t}; std(pupil_data_percent_signal_change)];
            
%             save([pupil_directory, 'pupil_data_percent_signal_change.mat'], 'pupil_data_percent_signal_change');
            pupil_data_percent_signal_change(:, 2) = PupilData(:, 7); % add triggers data
%             filename = [pupil_directory, 'pupil_data_percent_signal_change.mat'];
            setname = ['AB', num2str(p), '_', task{t}];
            EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',pupil_data_percent_signal_change','setname',setname,'srate',240);%,'pnts',0,'xmin',0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
            EEG = pop_chanevent(EEG, 2,'edge','leading','edgelen',0);
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = pop_epoch( EEG, {  '1'  }, [-1  6]); % epoch locked to cue onset
            EEG = pop_rmbase( EEG, [-200    0]); % subtract baseline pupil
            % delete epochs with artefacts as well as cue-only trials and error trials and trials after error trials
            EEG = pop_selectevent( EEG, 'omitepoch', RejectEpochsCorrect{t} ,'deleteevents','off','deleteepochs','on','invertepochs','off');
            
            % average pupil dilation from 1000 - 1500 ms after cue onset -
            % 240 sampling rate
            % pupil_avg_response{t} = squeeze(mean(EEG.data(1, 360:end, setdiff(1:60, RejectEpochsCorrect{t-1})), 2)); % time x trials 
            pupil_avg_response_tmp{t} = squeeze(mean(EEG.data(1, 481:600, :), 2)); % time x trials 
            pupil_response{t} = squeeze(EEG.data(1, :, :)); % time x trials
       
            
            % determine reaction time for each trial
            epoch_event=zeros(size(EEG.event, 2), 3);

            for e=1:size(EEG.event, 2)
                epoch_event(e, 1)=EEG.event(e).epoch;
                epoch_event(e, 2)=double(EEG.event(e).type);
                epoch_event(e, 3)=EEG.event(e).latency;
            end
            
            for r = 1:epoch_event(end, 1)
                events_per_epoch_index=[];
                events_per_epoch_index=find(epoch_event(:,1)==r);
                events_per_epoch=epoch_event(events_per_epoch_index, 2);
                latencies_per_epoch=epoch_event(events_per_epoch_index, 3);

                if events_per_epoch(2)==2 && events_per_epoch(3)==5 && length(events_per_epoch)==3
                    RT{t}=[RT{t}; (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/240]; 
%                     cue2target_interval = [cue2target_interval;(latencies_per_epoch(2)-latencies_per_epoch(1))*1000/240];
                end 
            end

        end
            % calculate std across trials, simple RT and gng task
            std_pupil_avg_response{grp, 1} = [std_pupil_avg_response{grp, 1}; std([pupil_avg_response_tmp{1};pupil_avg_response_tmp{2}])]; % simple RT
            std_pupil_avg_response{grp, 2} = [std_pupil_avg_response{grp, 2}; std([pupil_avg_response_tmp{3};pupil_avg_response_tmp{4}])]; % gng
            
            pupil_avg_response{grp, 1} = [pupil_avg_response{grp, 1}; mean([pupil_avg_response_tmp{1};pupil_avg_response_tmp{2}])]; % simple RT
            pupil_avg_response{grp, 2} = [pupil_avg_response{grp, 2}; mean([pupil_avg_response_tmp{3};pupil_avg_response_tmp{4}])]; % gng
            
            pupil_cuelocked_response{grp, 1} = [pupil_cuelocked_response{grp, 1}; mean([pupil_response{1}, pupil_response{2}], 2)'];
            pupil_cuelocked_response{grp, 2} = [pupil_cuelocked_response{grp, 2}; mean([pupil_response{3}, pupil_response{4}], 2)'];
            
            RT_single_trial{grp, 1}{count} = [RT{1}; RT{2}]; % group, task, participant
            RT_single_trial{grp, 2}{count} = [RT{3}; RT{4}];
            
            pupil_avg_single_trial{grp, 1}{count} = [pupil_avg_response_tmp{1}; pupil_avg_response_tmp{2}];
            pupil_avg_single_trial{grp, 2}{count} = [pupil_avg_response_tmp{3}; pupil_avg_response_tmp{4}];
            
            number_of_epochs{grp, 1}(count) = size([pupil_response{1}, pupil_response{2}], 2);
            number_of_epochs{grp, 2}(count) = size([pupil_response{3}, pupil_response{4}], 2);
            
            %% data to create graph showing the association between pupillary response and reaction time
%             pupil_response{t} = squeeze(EEG.data(1, :, :)); % time x trials

            pupil_single_trial_simpleRT = [pupil_response{1}, pupil_response{2}; RT_single_trial{grp, 1}{count}']';
            pupil_single_trial_gng = [pupil_response{3}, pupil_response{4}; RT_single_trial{grp, 2}{count}']';
            
            pupil_single_trial_simpleRT_sorted = sortrows(pupil_single_trial_simpleRT, size(pupil_single_trial_simpleRT, 2));
            pupil_single_trial_gng_sorted = sortrows(pupil_single_trial_gng, size(pupil_single_trial_gng, 2));
            
            data_simple_RT{grp}(count,1,:) = mean(pupil_single_trial_simpleRT_sorted(1:floor(size(pupil_single_trial_simpleRT_sorted, 1)/5), 1:end-1), 1);
            data_simple_RT{grp}(count,2,:) = mean(pupil_single_trial_simpleRT_sorted(floor(size(pupil_single_trial_simpleRT_sorted, 1)/5)+1:2*floor(size(pupil_single_trial_simpleRT_sorted, 1)/5), 1:end-1), 1);
            data_simple_RT{grp}(count,3,:) = mean(pupil_single_trial_simpleRT_sorted(2*floor(size(pupil_single_trial_simpleRT_sorted, 1)/5)+1:3*floor(size(pupil_single_trial_simpleRT_sorted, 1)/5), 1:end-1), 1);
            data_simple_RT{grp}(count,4,:) = mean(pupil_single_trial_simpleRT_sorted(3*floor(size(pupil_single_trial_simpleRT_sorted, 1)/5)+1:4*floor(size(pupil_single_trial_simpleRT_sorted, 1)/5), 1:end-1), 1);
            data_simple_RT{grp}(count,5,:) = mean(pupil_single_trial_simpleRT_sorted(4*floor(size(pupil_single_trial_simpleRT_sorted, 1)/5)+1:end, 1:end-1), 1);

            data_gng{grp}(count,1,:) = mean(pupil_single_trial_gng_sorted(1:floor(size(pupil_single_trial_gng_sorted, 1)/5), 1:end-1), 1);
            data_gng{grp}(count,2,:) = mean(pupil_single_trial_gng_sorted(floor(size(pupil_single_trial_gng_sorted, 1)/5)+1:2*floor(size(pupil_single_trial_gng_sorted, 1)/5), 1:end-1), 1);
            data_gng{grp}(count,3,:) = mean(pupil_single_trial_gng_sorted(2*floor(size(pupil_single_trial_gng_sorted, 1)/5)+1:3*floor(size(pupil_single_trial_gng_sorted, 1)/5), 1:end-1), 1);
            data_gng{grp}(count,4,:) = mean(pupil_single_trial_gng_sorted(3*floor(size(pupil_single_trial_gng_sorted, 1)/5)+1:4*floor(size(pupil_single_trial_gng_sorted, 1)/5), 1:end-1), 1);
            data_gng{grp}(count,5,:) = mean(pupil_single_trial_gng_sorted(4*floor(size(pupil_single_trial_gng_sorted, 1)/5)+1:end, 1:end-1), 1);
            
            
%         % graph of all trials per participant
%         data_correct = pupil_single_trial_gng(:, 1:end-1)'; x_axis=-.999:1/240:6;
%         for trl = 1:size(data_correct, 2)
%             clr(trl) = floor(255/size(data_correct, 2)*trl)/255;
%         end
%         color_index = randperm(size(data_correct, 2));
%         figure;
%         for trl = 1:size(data_correct, 2)
%     %         plot(x_axis, data_correct(16, :, trl), 'color', [clr(color_index(trl)) .5 .5]); hold on
%             plot(x_axis, data_correct(:, trl), 'color', [.5 clr(color_index(trl)) .5]); hold on
%         end
%         plot(x_axis, mean(data_correct(:, :), 2), 'k', 'LineWidth', 2); hold on 
%         plot(x_axis, zeros(length(x_axis), 1), '--k'); hold off
%         ax = gca; ax.FontSize = 20; ax.FontName = 'Arial'; ax.Color = 'none';
%         axis([0 4 -.3 .4])
%         xlabel('Time (s)', 'FontSize', 28, 'FontWeight','normal')
%         ylabel('Pupil response (%)', 'FontSize', 28, 'FontWeight','normal')
%         title(['AB', num2str(p), ' - Go-no/go'], 'FontSize', 28, 'FontWeight','normal')
  
    end
end

%% save variable
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability');
% save std_pupil_avg_response std_pupil_avg_response
% save std_pupil std_pupil
save pupil_cuelocked_response
% save pupil_avg_response pupil_avg_response
% save number_of_epochs number_of_epochs

%% plot data sorted by reaction time

plot_quintiles(data_simple_RT{1}, 'Simple RT - young');
plot_quintiles(data_gng{1}, 'Young');
plot_quintiles(data_simple_RT{2}, 'Simple RT - older');
plot_quintiles(data_gng{2}, 'Older');


%% within-subject correlation between pupil response amplitude and reaction time 
for grp = 1:2
    for task = 1:2
        for subj = 1:length(pupil_avg_single_trial{grp, task})
            [r,~,~,~,~,~] = skipped_correlation(pupil_avg_single_trial{grp, task}{subj},RT_single_trial{grp, task}{subj},0);
            coeff_robust{grp, task}(subj) = r.Pearson;
        end
    end
end

plot_all_data_2tasks(coeff_robust{1,1}, coeff_robust{1, 2},...
    coeff_robust{2, 1}, coeff_robust{2, 2}, 'Correlation\it r');

% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(coeff_robust{1, 2}, coeff_robust{2, 2}, 'Correlation\it r')


%% are coefficients different from zero
p_value = cell(2, 2); stats = cell(2, 2);
[h,p_value{1,1},ci,stats{1,1}] = ttest(coeff_robust{1,1})
[h,p_value{1,2},ci,stats{1,2}] = ttest(coeff_robust{1,2})
[h,p_value{2,1},ci,stats{2,1}] = ttest(coeff_robust{2,1})
[h,p_value{2,2},ci,stats{2,2}] = ttest(coeff_robust{2,2})

% are coefficients different across groups
p_value2 = cell(2, 1); stats2 = cell(2, 1);
[h,p_value2{1},ci,stats2{1}] = ttest2(coeff_robust{1,1}, coeff_robust{2,1}); % simple RT
[h,p_value2{2},ci,stats2{2}] = ttest2(coeff_robust{1,2}, coeff_robust{2,2}); % gng

% coefficients both groups together one-sample
p_value1 = cell(2, 1); stats1 = cell(2, 1);
[h,p_value1{1},ci,stats1{1}] = ttest([coeff_robust{1,1}, coeff_robust{2,1}]); % simple RT
[h,p_value1{2},ci,stats1{2}] = ttest([coeff_robust{1,2}, coeff_robust{2,2}]); % gng

%% creat excel file for SPSS analyses
% standard deviation over time, whole run
% standard deviation over trials - evoked cue-locked response
clear T_values
T_values(:, 1) = [ones(length(std_pupil_avg_response{1, 1}), 1); ones(length(std_pupil_avg_response{2, 1}), 1)*2];
for t = 1:2
    T_values(:, t+1) = [pupil_avg_response{1, t}; pupil_avg_response{2, t}];
    T_values(:, t+3) = [std_pupil_avg_response{1, t}; std_pupil_avg_response{2, t}];
end
clear std_pupil_total
std_pupil_total(:, 1) = [mean([std_pupil{1, 2}, std_pupil{1, 3}], 2); mean([std_pupil{2, 2}, std_pupil{2, 3}], 2)];
std_pupil_total(:, 2) = [mean([std_pupil{1, 4}, std_pupil{1, 5}], 2); mean([std_pupil{2, 4}, std_pupil{2, 5}], 2)];

T_values(:, 6:7) = std_pupil_total;

T_values(:, 8) = [coeff_robust{1,1}'; coeff_robust{2,1}'];
T_values(:, 9) = [coeff_robust{1,2}'; coeff_robust{2,2}'];


T = array2table(T_values, ...
    'VariableNames',{'Group' 'pupil_resp_simpleRT' 'pupil_resp_gng'  'pupil_resp_std_simpleRT' 'pupil_resp_std_gng' 'pupil_std_simpleRT' 'pupil_std_gng'...
    'corr_coef_simpleRT' 'corr_coef_gng'});

filename = 'pupil_std_500ms.xlsx';
writetable(T,filename,'Sheet',1,'Range','A1')

%% figures
% plot pupillary responses for both groups both tasks
% cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability');
% load pupil_cuelocked_response
% plot_data_mean_std_young_older(data_young, data_older, graph_title, y_axis_label)
plot_data_mean_std_young_older(pupil_cuelocked_response{1, 1}, pupil_cuelocked_response{2, 1}, 'Simple RT', 'Pupil response (%)')
plot_data_mean_std_young_older(pupil_cuelocked_response{1, 2}, pupil_cuelocked_response{2, 2}, '', 'Pupil response (%)')

%%
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability');
load std_pupil_avg_response; load pupil_avg_response;
% plot_all_data_2tasks(std_pupil_avg_response{1,1}, std_pupil_avg_response{1, 2},...
%     std_pupil_avg_response{2, 1}, std_pupil_avg_response{2, 2}, 'Pupil response SD');

% plot only go/no-go task
plot_all_data_onetask(std_pupil_avg_response{1, 2},  std_pupil_avg_response{2, 2}, 'Pupil response SD')

[H,P,CI,STATS] = ttest2(std_pupil_avg_response{1, 2},  std_pupil_avg_response{2, 2})

% load pupil_avg_response
% plot_all_data_2tasks(pupil_avg_response{1,1}, pupil_avg_response{1, 2},...
%     pupil_avg_response{2, 1}, pupil_avg_response{2, 2}, 'Pupil response (%)');

% plot only go/no-go task
plot_all_data_onetask(pupil_avg_response{1, 2}, pupil_avg_response{2, 2}, 'Pupil response (%)')
[H,P,CI,STATS] = ttest2(pupil_avg_response{1, 2}, pupil_avg_response{2, 2})

%% correlation between pupil and EEG spectral parameters - alpha power!
% participants used in eeg analyses
young_eeg = [4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older_eeg = [7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];

% ERP variability calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\eeglab_analysis_5_channellevel_AllChannels_CorrectNoErr_StDev.m
erp_variability_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability';
load([erp_variability_dir filesep 'ERP_avg_amp_stdev_young']); load([erp_variability_dir filesep 'ERP_avg_amp_stdev_older']); % task x participant x channel

% need to choose only the participants that are included in both analyses!
% simple RT
[R,P] = corrcoef(squeeze(ERP_avg_amp_stdev_young(1, find(young_eeg ~= setdiff(young_eeg, younger)), 16)), std_pupil_avg_response{1, 1})
figure; plot(squeeze(ERP_avg_amp_stdev_young(1, find(young_eeg ~= setdiff(young_eeg, younger)), 16)), std_pupil_avg_response{1, 1}, 'o');

[R,P] = corrcoef(squeeze(ERP_avg_amp_stdev_older(1, find(older_eeg ~= setdiff(older_eeg, older)), 16)), std_pupil_avg_response{2, 1})
figure; plot(squeeze(ERP_avg_amp_stdev_older(1, find(older_eeg ~= setdiff(older_eeg, older)), 16)), std_pupil_avg_response{2, 1}, 'o');

% gng
[R,P] = corrcoef(squeeze(ERP_avg_amp_stdev_young(2, find(young_eeg ~= setdiff(young_eeg, younger)), 16)), std_pupil_avg_response{1, 2})
figure; plot(squeeze(ERP_avg_amp_stdev_young(2, find(young_eeg ~= setdiff(young_eeg, younger)), 16)), std_pupil_avg_response{1, 2}, 'o');

[R,P] = corrcoef(squeeze(ERP_avg_amp_stdev_older(2, find(older_eeg ~= setdiff(older_eeg, older)), 16)), std_pupil_avg_response{2, 2})
figure; plot(squeeze(ERP_avg_amp_stdev_older(2, find(older_eeg ~= setdiff(older_eeg, older)), 16)), std_pupil_avg_response{2, 2}, 'o');

%% function to plot all data points 2 tasks 2 groups
function plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)

    % plot data for young group - simple RT and go/nogo task
    figure; box on; hold on
    
    % plot data for group 1
    yMean1=nanmean(data_grp1_task1); yMean2=nanmean(data_grp1_task2);
    y_se1 = nanstd(data_grp1_task1)/sqrt(length(data_grp1_task1)); y_se2 = nanstd(data_grp1_task2)/sqrt(length(data_grp1_task2));
    
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
    
    
    % group 2
        yMean1=nanmean(data_grp2_task1); yMean2=nanmean(data_grp2_task2);
    y_se1 = nanstd(data_grp2_task1)/sqrt(length(data_grp2_task1)); y_se2 = nanstd(data_grp2_task2)/sqrt(length(data_grp2_task2));
    
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
    
    % plot line on zero
    plot([0 6], zeros(2, 1), '--k')
    

    % axes('XColor','none');
    hold off;
    axis([0 6 -inf inf]);
    ax = gca;
    c = ax.Color;
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel=[1 2 1 2];
    xticks([1 2 4 5])
    ylabel(y_label_text, 'FontSize',32, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end

%% function plot data 2 groups 1 task
function plot_all_data_2groups(data_grp1, data_grp2, y_label_text, title_text)
       figure;
    % plot data for group 1
        yMean=nanmean(data_grp1);
    y_se = std(data_grp1)/sqrt(length(data_grp1));
    
        %plot the mean+-SEM box
        %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
%   Specify pos as a four-element vector of the form [x y w h] in data
%   units. The x and y elements determine the location and the w and h
%   elements determine the size. The function plots into the current axes
%   without clearing existing content from the axes.
    box on
    rectangle('Position',[1-0.3,yMean-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    hold on
    for y=1:length(data_grp1)
        plot(1+rand*0.2-0.1, data_grp1(y),'-o', 'color', [.8 .8 .8], ...
            'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end
    %plot line at zero    
    plot([0 4],[0 0],'--k','LineWidth', 1);
     %plot the mean line    
    plot([1-0.3 1+0.3],[yMean yMean] ,'Color','k','LineWidth',5);
   

    % group 2
    yMean=nanmean(data_grp2);
    y_se = std(data_grp2)/sqrt(length(data_grp2));
   rectangle('Position',[2-0.3,yMean-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1)
   
    for y=1:length(data_grp2)
        plot(2+rand*0.2-0.1, [data_grp2(y)] ,'-o', 'color', [1 .8 .8], ...
            'MarkerFaceColor',[1 .8 .8], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;  
    end

     %plot the mean line
    plot([2-0.3 2+0.3],[yMean yMean] ,'Color','k','LineWidth',5);
    

    % axes('XColor','none');
    hold off;
    axis([0 3 -inf inf]);
    ax = gca;
    c = ax.Color;
    ax.FontSize = 18;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel= [];
%     xticks([1 2])
    ylabel(y_label_text, 'FontSize', 24, 'FontWeight','normal')
    title(title_text, 'FontSize', 24, 'FontWeight','normal')
    
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

x_axis=-.999:1/240:6;
figure;
for x = 1.001:1/240:1.500 % grey background between 1000 amd 1500ms after ceu-onset
    plot([x x],[-.009 .099], 'color', [.8 .8 .8] ); hold on
end
hold on
% detection data
plot(x_axis, mean_data_young,'color', 'k', 'LineWidth',1.5)
hold on
jbfill(x_axis,mean_data_young+se_data_young, mean_data_young-se_data_young, 'k','k', 0.1)
% GNG data
hold on
plot(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
hold on
jbfill(x_axis, mean_data_older+se_data_older, mean_data_older-se_data_older,'r','r', 0.1)
hold on
plot(x_axis, zeros(length(x_axis), 1), '--k')
hold off
box off
ax = gca; ax.LineWidth = 2.5; 
% c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 24;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-.2 4 -.01 .1])
xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
ylabel(y_axis_label, 'FontSize', 32, 'FontWeight','normal')
title(graph_title, 'FontSize', 32, 'FontWeight','normal')

end


function plot_quintiles(data, title_text)

    Mean1= squeeze(mean(data(:,1,:), 1));
    Mean2= squeeze(mean(data(:,2,:), 1));
    Mean3= squeeze(mean(data(:,3,:), 1));
    Mean4= squeeze(mean(data(:,4,:), 1));
    Mean5= squeeze(mean(data(:,5,:), 1));

    SE1=(squeeze(std(data(:,1,:), 0, 1))/sqrt(size(data, 1)));
    SE2=(squeeze(std(data(:,2,:), 0, 1))/sqrt(size(data, 1)));
    SE3=(squeeze(std(data(:,3,:), 0, 1))/sqrt(size(data, 1)));
    SE4=(squeeze(std(data(:,4,:), 0, 1))/sqrt(size(data, 1)));
    SE5=(squeeze(std(data(:,5,:), 0, 1))/sqrt(size(data, 1)));

    colormap cool;%colormap summer;
    cmap = colormap;

    figure;
    xaxis = -.999:1/240:6;
    
%     for x = 1.001:1/240:1.500 % grey background between 1000 amd 1500ms after cue-onset
%         plot([x x],[-.0045 .0695], 'color', [.8 .8 .8] ); hold on
%     end
    
    % zero line
    plot(xaxis, zeros(length(xaxis)), 'k:');
    hold on
    plot( xaxis', Mean1, 'color', cmap(40,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean1+SE1)',(Mean1-SE1)', cmap(40,:), cmap(40,:), 1, 0.2)
    hold on
    plot( xaxis, Mean2, 'color', cmap(80,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean2+SE2)',(Mean2-SE2)', cmap(80,:), cmap(80,:), 1, 0.2)
    hold on
    plot( xaxis, Mean3, 'color', cmap(120,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean3+SE3)',(Mean3-SE3)', cmap(120,:), cmap(120,:), 1, 0.2)
    hold on
    plot( xaxis, Mean4, 'color', cmap(160,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean4+SE4)',(Mean4-SE4)', cmap(160,:), cmap(160,:), 1, 0.2)
    hold on
    plot( xaxis, Mean5, 'color', cmap(200,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean5+SE5)',(Mean5-SE5)', cmap(200,:), cmap(200,:), 1, 0.2)
    hold on
    % plot(xaxis, zeros(1, 840), ':k')
    hold off
    box off
    axis([0 2 -.005 .07]);
    % title('GNG', 'FontSize', 24, 'FontWeight','bold')
    % legend('RT1', 'RT2', 'RT3', 'RT4', 'RT5', 'location', 'northwest')
    ax = gca; ax.LineWidth = 2.5; 
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Pupil response (%)', 'FontSize', 32, 'FontWeight','normal')
    title(title_text, 'FontSize', 32, 'FontWeight','normal');
    
end



function plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)

    % plot data for young group - go/nogo task
    figure; box off; hold on
    
    % zero line
    plot([0 4], zeros(2, 1), 'k--');
    
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
    
    
    % group 2
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
    ax.LineWidth = 2.5; 
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