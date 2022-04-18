% study of the factors affecting the variability in the ERP amplitude at FCz in correct go trials
% factors studied: time-on-task; pre-stim spectral
% parameters (alpha power); pre-stim eeg phase
% amplitude
clear; close all;

% participants id
% eeg included participants
young_eeg=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older_eeg=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];

group = {young_eeg, older_eeg};
task={'D1', 'D2', 'G1', 'G2'};

% load epochs included in pre-stim spectral analysis - L:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\pre_stimulus_power_spectrum.m
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\epochs_prestim_older.mat');
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\epochs_prestim_young.mat');

% alpha power in POz calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\foof_single_trial_spectrum_POz.m
alpha_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\single_trial_pre_stim_spectra_POz';

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

betas_young = cell(2, 1); betas_older = cell(2, 1); erp_resd_alpha_bslnphase_std = cell(2, 2); eeg_model_rsquared = cell(2, 2);
for grp = 1:2
    for p = group{grp}

 
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        
        % spectral properties in pre-cue time window (baseline)-alpha
        % power and exponent calculated in 
        % G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\foof_single_trial_spectrum_POz.m
        epochs_prestim_spectra = {}; alpha_power = {};% beta_power = {}; exponent_prestim = {}; 
        if grp == 1
%             addpath([alpha_dir filesep 'young']);
            load([alpha_dir filesep 'young' filesep 'part_id_' num2str(p)]);
            for t = 1:4
                epochs_prestim_spectra{t} = epochs_prestim_young{young_eeg == p, t+1}; % excluding passive task - first run
                alpha_power{t} = oscillations_power{t+1}(:, 2);
%                 beta_power{t} = oscillations_power{t+1}(:, 3);
%                 exponent_prestim{t} = exponent{t+1};
            end
        else
%             addpath([alpha_dir filesep 'older']);
            load([alpha_dir filesep 'older' filesep 'part_id_' num2str(p)]);
            for t = 1:4
                epochs_prestim_spectra{t} = epochs_prestim_older{older_eeg == p, t+1};
                alpha_power{t} = oscillations_power{t+1}(:, 2);
%                 beta_power{t} = oscillations_power{t+1}(:, 3);
%                 exponent_prestim{t} = exponent{t+1};
            end
        end
     
%         %% erp data - dependent variable
        eeg_directory = strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\EEG\');
% 
        % create variable with data electrodes X time frames X trials -
        % starting at cue onset
        epochs_erp = {}; FCz_erp = {}; 
       for t = 1:4
           % load data created in G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\eeglab_analysis_3.m
           filename=strcat('AB', num2str(p), '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej_epochs_baselined.set');
           EEG = pop_loadset(filename, eeg_directory);
           
           % determine original epoch number to match with pre-stim alpha
           % for each original epoch find the event original event (cue) number 
            a=0; trial_numbers=[]; trial_numbers(:, 1)=1:60; % there are 60 trials in each run - original trial number
            for i=1:size(EEG.urevent, 2)
                if EEG.urevent(i).type == 1
                    a=a+1;
                    trial_numbers(a, 2) = i; % event number on original file without periods with artefacts removed
                elseif strcmp(EEG.urevent(i).type, {'boundary'})
                    trial_numbers_boundary = i;
                end
            end
           epochs_erp{t} = zeros(size(EEG.epoch, 2), 2); % epoch number x epoch original
           epochs_erp{t}(:, 1) = 1:size(EEG.epoch, 2);
           for epoch = 1:size(EEG.epoch, 2)
                for i = 1:size(EEG.event, 2)
                    if EEG.event(i).epoch == epoch && strcmp(EEG.event(i).type, '1')
                        epochs_erp{t}(epoch, 2) = trial_numbers((trial_numbers(:, 2) == EEG.event(i).urevent), 1);
                    end
                end
           end

           % FCz ERP data
           for chan = 1:length(EEG.chanlocs)
               if strcmp(EEG.chanlocs(chan).labels, 'FCz')
                    FCz_erp{t} = squeeze(EEG.data(chan, :, :)); % time X trials
               end
           end
       end
       
       
       % ERP baseline - phase and amplitude envelope of slow fluctuations
       for t = 1:4
             % clear eeglab
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            %load file - from eeglab_analysis_2b_synchronize_pupil
            % _RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica_ArtICsRem_Filt0_1_35Hz_Pupil.set
            % manual/visual removal of periods with artifacts
            % save as ..._Filt0_1_35Hz_Pupil_ManualArtRej.set
            if p == 77 && ismember(t, [1, 2])
                filename=strcat('AB', num2str(p), '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej_HEO.set');
            else
                filename=strcat('AB', num2str(p), '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej.set');
            end
            EEG = pop_loadset(filename, eeg_directory);
            % find which chanel is FCz
            for chan = 1:length(EEG.chanlocs)
               if strcmp(EEG.chanlocs(chan).labels, 'FCz')
                    c = chan;
               end
           end
            
            % filter data low pass 2 Hz
            EEG = pop_eegfiltnew(EEG, 'hicutoff',2,'plotfreqz',0);
            % hilbert transform of FCz
            signal_phase = [];
            y = hilbert(EEG.data(c, :));
            signal_phase = angle(y); % extract instantaneous phase
            signal_amplitude = abs(y); % extract instantaneous amplitude envelope

            EEG.data(c, :) = signal_phase;
            EEG.data(c+1, :) = signal_amplitude;
            % epoch data - 3 sec before cue onset
            EEG = pop_epoch( EEG, {  '1'  }, [-3.5  0]);
            % FCz ERP data

           % phase and amplitude at cue-onset
           FCz_erp_bsln{t} = squeeze(EEG.data(c, end, :));
           FCz_erp_bslnamp{t} = squeeze(EEG.data(c+1, end, :));
       end

       %% match epochs of alpha, and erp
%        epochs_prestim_spectra
%        epochs_erp
       included_epochs = {};
       for t = 1:4
            included_epochs{t} = mintersect(epochs_prestim_spectra{t}(:, 2), epochs_erp{t}(:, 2));

            % delete epochs not to include
            % alpha, beta and exponent
            index2delete = []; epochs2delete = setdiff(epochs_prestim_spectra{t}(:, 2), included_epochs{t});
            if ~isempty(epochs2delete)
                for epoch = 1:length(epochs2delete)
                    index2delete = [index2delete; find(epochs_prestim_spectra{t}(:, 2) == epochs2delete(epoch))];
                end
                alpha_power{t}(index2delete ) = []; 
%                 beta_power{t}(index2delete ) = []; 
%                 exponent_prestim{t}(index2delete ) = []; 
                FCz_erp_bsln{t}(index2delete ) = []; % calculated in the same epochs as alpha
                FCz_erp_bslnamp{t}(index2delete ) = []; % calculated in the same epochs as alpha
            end
            % erp
            index2delete = []; epochs2delete = setdiff(epochs_erp{t}(:, 2), included_epochs{t});
            if ~isempty(epochs2delete)
                for epoch = 1:length(epochs2delete)
                    index2delete = [index2delete; find(epochs_erp{t}(:, 2) == epochs2delete(epoch))];
                end
                FCz_erp{t}(:, index2delete) = []; 
            end
       end
       

         %% create variables for model
         task_runs = [1 2; 3 4];
         for tsk = 1:2 % simple RT and gng
             correct_go_trials_data = [];  categorical_variable = []; 
             % dependent variable FCz ERP from 1000:1500 ms from cue onset
            correct_go_trials_data(:, 1)=[mean(FCz_erp{task_runs(tsk, 1)}(601:850, :), 1)'; mean(FCz_erp{task_runs(tsk, 2)}(601:850, :), 1)'];
            % categorical variable with run 
            categorical_variable(:,1)=[zeros(size(FCz_erp{task_runs(tsk, 1)}, 2), 1); ones(size(FCz_erp{task_runs(tsk, 2)}, 2), 1)];

        % create continuous variable with time-on-task; baseline pupil ;alpha baseline - total 3 continuous variables
            total_length = length(included_epochs{task_runs(tsk, 1)}) + length(included_epochs{task_runs(tsk, 2)});
            continuous_variable = zeros(total_length, 3);
            % time-on-task
            % z-score time on task variable
    %         continuous_variable(: ,1) = zscore([included_epochs{1}; included_epochs{2}; included_epochs{3}; included_epochs{4}]);
            continuous_variable(: ,1) = [included_epochs{task_runs(tsk, 1)}; included_epochs{task_runs(tsk, 2)}];

            % pre-stimulus alpha power
    %         continuous_variable(:, 2) = zscore([alpha_power{1}; alpha_power{2}; alpha_power{3}; alpha_power{4}]);
            continuous_variable(:, 2) = [alpha_power{task_runs(tsk, 1)}; alpha_power{task_runs(tsk, 2)}];

             % erp baseline phase - cos(angle)
            continuous_variable(:, 3) = [FCz_erp_bslnamp{task_runs(tsk, 1)}.*cos(FCz_erp_bsln{task_runs(tsk, 1)}); FCz_erp_bslnamp{task_runs(tsk, 2)}.*cos(FCz_erp_bsln{task_runs(tsk, 2)})];

            g = fitglm([categorical_variable, continuous_variable], correct_go_trials_data, 'Categorical',1);
            eeg_model_rsquared{tsk, grp} = [eeg_model_rsquared{tsk, grp}; g.Rsquared.Ordinary, g.Rsquared.Adjusted];

            %% save coefficients
            if grp == 1
                betas_young{tsk} = [betas_young{tsk}; g.Coefficients.Estimate'];
            else
                betas_older{tsk} = [betas_older{tsk}; g.Coefficients.Estimate'];
            end

            % residuals after regressing all factors out
            erp_resd_alpha_bslnphase_std{tsk, grp} = [erp_resd_alpha_bslnphase_std{tsk, grp}; std(g.Residuals.Raw)];
        end 
    end
end
save eeg_modelbslnamp_rsquared_per_task eeg_model_rsquared
save betas_time_alphaPOz_EEGbslnamp_per_task_young betas_young
save betas_time_alphaPOz_EEGbslnamp_per_task_older betas_older

%% stats
% cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability')
% load betas_time_alphaPOz_EEGbslnampphase_per_task_young; load betas_time_alphaPOz_EEGbslnampphase_per_task_older;
% % task, run, time-on-task, alpha, pupil, eeg bsln
% [h,p,ci,stats] = ttest(betas_young)
% 
% [h,p,ci,stats] = ttest(betas_older)
% 
% [h,p,ci,stats] = ttest2(betas_young, betas_older)
% 
% [h,p,ci,stats] = ttest([betas_young; betas_older])


%% model r squared 
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability')
load eeg_modelbslnamp_rsquared_per_task % task x group
% eeg_model_rsquared{tsk, grp}
plot_all_data_2groups(eeg_model_rsquared{2, 1}(:, 1), eeg_model_rsquared{2, 2}(:,1),'Model R squared', '')

% comparison of model r squared with and without sin(theta)
rsquared_withoutsintheta = eeg_model_rsquared;

load eeg_model_rsquared_per_task; % full model calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\CueLockFCzvsTimeOnTask_AlphaPOz_FCzbslnampphase_per_task.m

[h,p,ci,stats] = ttest([rsquared_withoutsintheta{2, 1}(:, 1); rsquared_withoutsintheta{2, 2}(:, 1)],...
    [eeg_model_rsquared{2, 1}(:, 1); eeg_model_rsquared{2, 2}(:,1)]) 

mean(rsquared_withoutsintheta{2, 1}(:, 1))
std(rsquared_withoutsintheta{2, 1}(:, 1))
mean(rsquared_withoutsintheta{2, 2}(:, 1))
std(rsquared_withoutsintheta{2, 2}(:, 1))


plot_all_data_2tasks(eeg_model_rsquared{2, 1}(:, 1), rsquared_withoutsintheta{2, 1}(:, 1), eeg_model_rsquared{2, 2}(:, 1), rsquared_withoutsintheta{2, 2}(:, 1), 'Model \itR\rm2')

%% erp_residuals_total_std{grp}
plot_all_data_2groups(erp_resd_alpha_bslnphase_std{1}, erp_resd_alpha_bslnphase_std{2}, 'ERP residuals std', '')

% % % create excel file for SPSS analysis
% % clear T_values
% % T_values(:, 1) = [ones(size(betas_young, 1), 1); ones(size(betas_older, 1), 1)*2];
% % T_values(:, 2:8) =  [betas_young(:, 2:end); betas_older(:, 2:end)];
% % 
% % all_column_names = {'group','task','run','Time_on_task' 'AlphaPower' 'Prestim_pupil' 'Prestim_EEGcos' 'Prestim_EEGsin'};
% % 
% % T = array2table(T_values, ...
% %     'VariableNames',all_column_names);
% % 
% % filename = 'betas_time_alphaPOz_pupil_eegbslnphase.xlsx';
% % writetable(T,filename,'Sheet',1,'Range','A1')


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