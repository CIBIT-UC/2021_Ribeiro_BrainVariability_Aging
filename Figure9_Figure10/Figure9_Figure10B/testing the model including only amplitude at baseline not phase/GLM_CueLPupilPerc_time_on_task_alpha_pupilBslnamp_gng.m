% GLMs to investigate effect of run, time-on-run, pre-stim alpha and
% pre-stim pupil ou cue-locked pupil response
clear; close all;
% participants id
% eeg included participants - for alpha analysis
young_eeg=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older_eeg=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];

% pupil included participants
young_pupil = [4	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46 50 51 53 54 56	59	62 66	68	72	74	76  78  80 81   82  84  85];
older_pupil = [7  8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77 79   83  86];

% find participants with eeg and pupil data
young = mintersect(young_eeg, young_pupil);
older = mintersect(older_eeg, older_pupil);

group = {young, older};
task={'D1', 'D2', 'G1', 'G2'};

% load epochs included in pre-stim spectral analysis - L:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\pre_stimulus_power_spectrum.m
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\epochs_prestim_older.mat');
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\epochs_prestim_young.mat');

% alpha power in POz calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\foof_single_trial_spectrum_POz.m
alpha_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\single_trial_pre_stim_spectra_POz';

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

betas_young = []; betas_older = []; %R_alpha_pupil_young = []; R_alpha_pupil_older = [];
alpha_std_young = []; alpha_std_older = []; pupil_std_young = []; pupil_std_older = []; 
for grp = 1:2
    part = 0;
    for p = group{grp}
        part = part + 1;
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        % eeg directory
        eeg_dir = strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\EEG\');
        % pupil rejected epochs - these were manually/visually rejected
        pupil_directory=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\');
        % load rejected epochs needs data from:
        % G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\pupil_analysis_3_epoching.m
        %  G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\pupil_analysis_3a_behavior_per_run_b4artifacts.m
        % G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\pupil_analysis_3b_RejectedEpochs_Go_CueOnly_SaveRunData.m
        load([pupil_directory, 'RejectEpochsCorrectGoCueOnly']);
            pupil_baseline = {}; pupil_data = {}; time_on_task = {};
        for t = 3:4 % task runs excluding first run = passive task
                % calculate pupil go trials data
                filename=strcat(pupil_directory, 'AB', num2str(p), '_', task(t), '_PupilData.mat');
                load(filename{1})
                % use smoothed pupil diameter data with blink artefacts interpolated
                % measured in horizontal direction (x)
                pupil_data_percent_signal_change = PupilData(:, 3)/mean(PupilData(:, 3));
                pupil_data_percent_signal_change(:, 2) = PupilData(:, 7); % add triggers data
                save([pupil_directory, 'pupil_data_percent_signal_change.mat'], 'pupil_data_percent_signal_change');
                filename = [pupil_directory, 'pupil_data_percent_signal_change.mat'];
                setname = ['AB', num2str(p), '_', task{t}];
                EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',filename,'setname',setname,'srate',240);%,'pnts',0,'xmin',0);
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
                EEG = pop_chanevent(EEG, 2,'edge','leading','edgelen',0);
                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = pop_epoch( EEG, {  '1'  }, [-.2  1.5]); % epoch locked to cue onset
                EEG = pop_rmbase( EEG, [-200    0]); % subtract baseline pupil
                % average across time from 1000 to 1500 ms after cue onset
                % delete epochs with artefacts as well as cue-only trials and error trials and trials after error trials
                pupil_data{t} = squeeze(mean(EEG.data(1, 289:end, setdiff(1:60, RejectEpochsCorrectGoCueOnly{t})), 2)); % channel x time x trials

    %           % create variable with trial number for each run - time on task
    %           = trial number
                time_on_task{t}(:,1)=setdiff(1:60, RejectEpochsCorrectGoCueOnly{t})';% exclude trials that had artefacts
        end

        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        % calculate phase of pupil fluctuations at cue-onset
        for t = 3:4 % task runs excluding first run = passive task
                % calculate pupil go trials data
                filename=strcat(pupil_directory, 'AB', num2str(p), '_', task(t), '_PupilData.mat');
                load(filename{1})
                % use smoothed pupil diameter data with blink artefacts interpolated
                % measured in horizontal direction (x)
                pupil_data_percent_signal_change = PupilData(:, 3)/mean(PupilData(:, 3));
                pupil_data_percent_signal_change(:, 2) = PupilData(:, 7); % add triggers data
                save([pupil_directory, 'pupil_data_percent_signal_change.mat'], 'pupil_data_percent_signal_change');
                filename = [pupil_directory, 'pupil_data_percent_signal_change.mat'];
                setname = ['AB', num2str(p), '_', task{t}];
                EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',filename,'setname',setname,'srate',240);%,'pnts',0,'xmin',0);
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
                EEG = pop_chanevent(EEG, 2,'edge','leading','edgelen',0);
                % filter data between .1 and .9 Hz
                EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'hicutoff',0.9);
                % hilbert transform
                y = hilbert(EEG.data);
                signal_phase = angle(y); % extract instantaneous phase
                signal_amplitude = abs(y); % extract amplitude envelop
                EEG.data(1, :) = signal_phase;
                EEG.data(2, :) = signal_amplitude;
                EEG = pop_epoch( EEG, {  '1'  }, [-.2  1.5]); % epoch locked to cue onset
                pupil_bsln_phase{t} = squeeze(EEG.data(1, 48, setdiff(1:60, RejectEpochsCorrectGoCueOnly{t}))); % phase at cue-onset for each trial
                pupil_bsln_amp{t} = squeeze(EEG.data(2, 48, setdiff(1:60, RejectEpochsCorrectGoCueOnly{t}))); % amplitude envelop at cue-onset for each trial
        end

        %% spectral properties in pre-cue time window (baseline)-alpha
        % power and exponent calculated in 
        % G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\foof_single_trial_spectrum_POz.m
        epochs_prestim_spectra = {}; alpha_power = {};% beta_power = {}; exponent_prestim = {}; 
        if grp == 1
%             addpath([alpha_dir filesep 'young']);
            load([alpha_dir filesep 'young' filesep 'part_id_' num2str(p)]);
            for t = 3:4
                epochs_prestim_spectra{t} = epochs_prestim_young{young_eeg == p, t+1}; % excluding passive task - first run
                alpha_power{t} = oscillations_power{t+1}(:, 2);
%                 beta_power{t} = oscillations_power{t+1}(:, 3);
%                 exponent_prestim{t} = exponent{t+1};
            end
        else
%             addpath([alpha_dir filesep 'older']);
            load([alpha_dir filesep 'older' filesep 'part_id_' num2str(p)]);
            for t = 3:4
                epochs_prestim_spectra{t} = epochs_prestim_older{older_eeg == p, t+1};
                alpha_power{t} = oscillations_power{t+1}(:, 2);
%                 beta_power{t} = oscillations_power{t+1}(:, 3);
%                 exponent_prestim{t} = exponent{t+1};
            end
        end  
        
           
       %% match epochs of pupil, alpha
%        pupil_epochs
%        epochs_prestim_spectra
%        epochs_erp
       included_epochs = {}; pupil_epochs = time_on_task;
       for t = 3:4
            included_epochs{t} = mintersect(pupil_epochs{t}, epochs_prestim_spectra{t}(:, 2));

            % delete epochs not to include
            % pupil
            index2delete = []; epochs2delete = setdiff(pupil_epochs{t}, included_epochs{t});
            if ~isempty(epochs2delete)
                for epoch = 1:length(epochs2delete)
                    index2delete = [index2delete; find(pupil_epochs{t} == epochs2delete(epoch))];
                end
                pupil_bsln_phase{t}(index2delete) = []; 
                pupil_bsln_amp{t}(index2delete) = []; 
                pupil_data{t}(index2delete) = [];
                time_on_task{t}(index2delete) = [];
            end
            % alpha
            index2delete = []; epochs2delete = setdiff(epochs_prestim_spectra{t}(:, 2), included_epochs{t});
            if ~isempty(epochs2delete)
                for epoch = 1:length(epochs2delete)
                    index2delete = [index2delete; find(epochs_prestim_spectra{t}(:, 2) == epochs2delete(epoch))];
                end
                alpha_power{t}(index2delete ) = []; 
            end
       end
       
            %% create variables for model
     correct_go_trials_data = [];  categorical_variable = []; 
     % dependent variable pupil response from 1000:1500 ms from cue onset
        correct_go_trials_data(:, 1)=[pupil_data{1}; pupil_data{2}; pupil_data{3}; pupil_data{4}];
        % categorical variable with run
        categorical_variable(:,1)=[zeros(length(pupil_data{1}), 1); ones(length(pupil_data{2}), 1);...
         zeros(length(pupil_data{3}), 1); ones(length(pupil_data{4}), 1)];

    % create continuous variable with time-on-task; alpha baseline; pupil phase cos;
        total_length = length(included_epochs{1}) + length(included_epochs{2}) + length(included_epochs{3}) + length(included_epochs{4});
        continuous_variable = zeros(total_length, 3);
        % time-on-task
        % z-score time on task variable
%         continuous_variable(: ,1) = zscore([included_epochs{1}; included_epochs{2}; included_epochs{3}; included_epochs{4}]);
        % WITHOUT ZSCORE
        continuous_variable(: ,1) = [included_epochs{1}; included_epochs{2}; included_epochs{3}; included_epochs{4}];

        % pre-stimulus alpha power
%         continuous_variable(:, 2) = zscore([alpha_power{1}; alpha_power{2}; alpha_power{3}; alpha_power{4}]);
        continuous_variable(:, 2) = [alpha_power{1}; alpha_power{2}; alpha_power{3}; alpha_power{4}];

         % pupil baseline - phase of slow fluctuatons - cos sin
        continuous_variable(:, 3) = [pupil_bsln_amp{1}.*cos(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*cos(pupil_bsln_phase{2}); pupil_bsln_amp{3}.*cos(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*cos(pupil_bsln_phase{4})];
%         continuous_variable(:, 4) = [pupil_bsln_amp{1}.*sin(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*sin(pupil_bsln_phase{2}); pupil_bsln_amp{3}.*sin(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*sin(pupil_bsln_phase{4})];
       
        g = fitglm([categorical_variable, continuous_variable], correct_go_trials_data, 'Categorical',1);
        pupil_model_rsquared{grp}(part, :) = [g.Rsquared.Ordinary, g.Rsquared.Adjusted];
        pupil_resid_std{grp}(part) = std(g.Residuals.Raw);
        
        
        
        % save coefficients
        if grp == 1
            betas_young = [betas_young; g.Coefficients.Estimate'];
        else
            betas_older = [betas_older; g.Coefficients.Estimate'];

        end
        
    end
end

save betas_alpha_pupilamp_gng_young betas_young
save betas_alpha_pupilamp_gng_older betas_older

save pupil_modelbslnamp_rsquared_gng pupil_model_rsquared

%% stats
% task, run, time-on-task, alpha bsln, pupil bsln, eeg bsln


[h,p,ci,stats] = ttest(betas_young)

[h,p,ci,stats] = ttest(betas_older)


[h,p,ci,stats] = ttest2(betas_young, betas_older)

[h,p,ci,stats] = ttest([betas_young; betas_older])



[h,p,ci,stats] = ttest2(alpha_std_young, alpha_std_older)

[h,p,ci,stats] = ttest2(pupil_std_young, pupil_std_older)

%% graphs of coefficients

y_label_text = 'Coefficients';
title_text = {'Run' 'Time-on-task' 'Alpha power' 'Cos(\theta)' 'Sin(\theta)'}; 
    
for coef = 1:5
    % plot_all_data_2groups(data_grp1, data_grp2, y_label_text)
    plot_all_data_2groups(betas_young(: , coef+1), betas_older(: , coef+1), y_label_text, title_text(coef))
end

%% plot model Rsquared
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability')
load pupil_modelbslnamp_rsquared_gng
% plot_all_data_2groups(data_grp1, data_grp2, y_label_text)
plot_all_data_2groups(pupil_model_rsquared{1}(:, 1), pupil_model_rsquared{2}(:,1),'Model R squared', '')

plot_all_data_2groups(pupil_model_rsquared{1}(:, 2), pupil_model_rsquared{2}(:, 1),'Model R squared adjusted', '')

[h,p,ci,stats] = ttest2(pupil_model_rsquared{1}(:, 1), pupil_model_rsquared{2}(:,1))

[h,p,ci,stats] = ttest2(pupil_model_rsquared{1}(:, 2), pupil_model_rsquared{2}(:,2))


% comparison of model r squared with and without sin(theta)
rsquared_withoutsintheta = pupil_model_rsquared;

load pupil_model_rsquared_gng; % full model calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\GLM_CueLPupilPerc_time_on_task_alpha_pupilBslnphase_gng.m

[h,p,ci,stats] = ttest([rsquared_withoutsintheta{1}(:, 1); rsquared_withoutsintheta{2}(:, 1)],...
    [pupil_model_rsquared{1}(:, 1); pupil_model_rsquared{2}(:,1)])

mean(rsquared_withoutsintheta{1}(:, 1))
std(rsquared_withoutsintheta{1}(:, 1))
mean(rsquared_withoutsintheta{2}(:, 1))
std(rsquared_withoutsintheta{2}(:, 1))

%~ plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(pupil_model_rsquared{1}(:, 1), rsquared_withoutsintheta{1}(:, 1),...
    pupil_model_rsquared{2}(:,1), rsquared_withoutsintheta{2}(:, 1),  'Model \itR\rm2')

%% graph of correlation across independent variables
% plot_all_data_2groups(R_alpha_pupil_young(:, 1), R_alpha_pupil_older(:, 1), 'Correlation r','Alpha - time-on-task')
% 
% plot_all_data_2groups(R_alpha_pupil_young(:, 2), R_alpha_pupil_older(:, 2), 'Correlation r', 'Pupil - time-on-task')

%% alpha and pupil std
plot_all_data_2groups(alpha_std_young, alpha_std_older, 'Pre-stim alpha std', '')

plot_all_data_2groups(pupil_std_young, pupil_std_older, 'Pre-stim pupil std', '')

% %% create excel file for SPSS analysis
% clear T_values
% T_values(:, 1) = [ones(size(betas_young, 1), 1); ones(size(betas_older, 1), 1)*2];
% T_values(:, 2:6) =  [betas_young(:, 2:end); betas_older(:, 2:end)];
% 
% all_column_names = {'group','task','run','Time_on_task' 'AlphaPower' 'Prestim_pupil' 'alpha_timeontask' 'pupil_timeontask'};
% 
% T = array2table(T_values, ...
%     'VariableNames',all_column_names);
% 
% filename = 'betas_time_alphaPOz_pupil.xlsx';
% writetable(T,filename,'Sheet',1,'Range','A1')


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
%     %plot line at zero    
%     plot([0 4],[0 0],'--k','LineWidth', 1);
     %plot the mean line    
    plot([1-0.3 1+0.3],[yMean yMean] ,'Color','k','LineWidth',5);
   

    %% group 2
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
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
%     ax.XTickLabel= [];
    ax.XAxis.FontSize = 32;
    ax.XTickLabel= [{'' 'Young' 'Older' ''}];
%     xticks([1 2])
    ylabel(y_label_text, 'FontSize', 32, 'FontWeight','normal')
    title(title_text, 'FontSize', 32, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end

% function to plot all data points 2 tasks 2 groups
function plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)

    % plot data for young group - simple RT and go/nogo task
    figure; box off; hold on
    
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
        plot([4 5]+rand*0.2-0.1, [data_grp2_task1(y) data_grp2_task2(y)] ,'-o', 'color', [1 .5 .5], ...
            'MarkerFaceColor',[1 .5 .5], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1.5);
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
    ax.LineWidth = 2.5;
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel=[];
    xticks([1 2 4 5])
    ylabel(y_label_text, 'FontSize',32, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end