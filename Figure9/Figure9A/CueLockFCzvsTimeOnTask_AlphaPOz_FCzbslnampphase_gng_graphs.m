% study of the factors affecting the variability in the ERP amplitude at FCz in correct go trials
% factors studied: time-on-task; pre-stim pupil
% diameter; pre-stim spectral parameters (alpha power); pre-stim eeg
% amplitude
clear; close all;

% participants id
% eeg included participants
young_eeg=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older_eeg=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];

% pupil included participants
young_pupil = [4	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46 50 51 53 54 56	59	62 66	68	72	74	76  78  80 81   82  84  85];
older_pupil = [7  8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77 79   83  86];

% % heart rate included participants
% young_heart = setdiff([4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85], [10 17 19 35 61]);
% older_heart=setdiff([7	8	11	12	14	17	19	20	21	22	23	32	35	37	38  41	43	47	48	49	52	55	57	58	60	61	63	64	65	67	69	70	71	73	75	77	79	83  86], [10 17 19 35 61]);

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

% % dir with heart data
% heart_directory='G:\ProjectAgingNeuromodulation\AuditoryResearch\HeartRateAnalysis\TaskRelatedModulation\QRS_files_WithoutCorrection4EctopicBeats_mat';

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

betas_young = []; betas_older = [];
for grp = 1:2
    count = 0;
    for p = group{grp}
        count = count + 1;

 
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        
        % spectral properties in pre-cue time window (baseline)-alpha
        % power and exponent calculated in 
        % G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\foof_single_trial_spectrum_POz.m
        epochs_prestim_spectra = {}; alpha_power = {};% beta_power = {}; exponent_prestim = {}; 
        if grp == 1
%             addpath([alpha_dir filesep 'young']);
            load([alpha_dir filesep 'young' filesep 'part_id_' num2str(p)]);
            for t = 3:4%1:4
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
     
%         %% erp data - dependent variable
        eeg_directory = strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\EEG\');
% 
        % create variable with data electrodes X time frames X trials -
        % starting at cue onset
        epochs_erp = {}; FCz_erp = {}; 
       for t = 3:4
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
       
       
       % ERP baseline - phase of slow fluctuations
       for t = 3:4
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
            
%             f = figure;
%             time = .002:.002:60;
%             plot(time, EEG.data(c, 60*500+1:120*500), 'LineWidth', 1.5, 'color', [0, 0.4470, 0.7410]); hold on
%             plot(time, zeros(1, length(time)), '--k'); hold on
%             plot(time, signal_phase(60*500+1:120*500)*3,'r', 'LineWidth', 1)%, 'color', 	[0.8500, 0.3250, 0.0980])
%             ax = gca;
%             ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
%             xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
%             ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
%             f.Position = [10 10 1300 400]; 
            

            EEG.data(c, :) = signal_phase;
            EEG.data(c+1, :) = signal_amplitude;
            % epoch data - 3 sec before cue onset
            EEG = pop_epoch( EEG, {  '1'  }, [-3.5  0]);
            % FCz ERP data

           % phase and amplitude at cue-onset
           FCz_erp_bsln{t} = squeeze(EEG.data(c, end, :));
           FCz_erp_bslnamp{t} = squeeze(EEG.data(c+1, end, :));
       end
       
       
       
       % ERP baseline - baseline amplitude after filtering
       for t = 3:4
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
            % filter data low pass 2 Hz
            EEG = pop_eegfiltnew(EEG, 'hicutoff',2,'plotfreqz',0);
            % epoch data - 3 sec before cue onset
            EEG = pop_epoch( EEG, {  '1'  }, [-3.5  0]);
            % FCz ERP data
           for chan = 1:length(EEG.chanlocs)
               if strcmp(EEG.chanlocs(chan).labels, 'FCz')
                    FCz_erp_bsln_tmp = squeeze(EEG.data(chan, :, :)); % time X trials
               end
           end
           % amplitude at cue-onset
           FCz_erp_bsln_amplitude{t} = FCz_erp_bsln_tmp(end, :);
       end


       %% match epochs of alpha, and erp
%        epochs_prestim_spectra
%        epochs_erp
       included_epochs = {};
       for t = 3:4
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
                FCz_erp_bsln_amplitude{t}(index2delete ) = [];
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
       

     % create variables for model
     correct_go_trials_data = [];  categorical_variable = []; 
     % FCz ERP = time x trials
        correct_go_trials_data = [FCz_erp{1}'; FCz_erp{2}'; FCz_erp{3}'; FCz_erp{4}'];
        % categorical variable with task and run - one in each column
        categorical_variable(:,1)=[zeros(size(FCz_erp{1}, 2)+size(FCz_erp{2}, 2), 1); ones(size(FCz_erp{3}, 2)+size(FCz_erp{4}, 2), 1)];
        categorical_variable(:,2)=[zeros(size(FCz_erp{1}, 2), 1); ones(size(FCz_erp{2}, 2), 1);...
         zeros(size(FCz_erp{3}, 2), 1); ones(size(FCz_erp{4}, 2), 1)];

    % create continuous variable with time-on-task; baseline pupil ;alpha baseline - total 3 continuous variables
        total_length = length(included_epochs{1}) + length(included_epochs{2}) + length(included_epochs{3}) + length(included_epochs{4});
        continuous_variable = zeros(total_length, 3);
        % time-on-task
        % z-score time on task variable
%         continuous_variable(: ,1) = zscore([included_epochs{1}; included_epochs{2}; included_epochs{3}; included_epochs{4}]);
        continuous_variable(: ,1) = [included_epochs{1}; included_epochs{2}; included_epochs{3}; included_epochs{4}];

        % pre-stimulus alpha power
%         continuous_variable(:, 2) = zscore([alpha_power{1}; alpha_power{2}; alpha_power{3}; alpha_power{4}]);
        continuous_variable(:, 2) = [alpha_power{1}; alpha_power{2}; alpha_power{3}; alpha_power{4}];
        
         % erp baseline phase - cos(angle)
        continuous_variable(:, 3) = [FCz_erp_bslnamp{1}.*cos(FCz_erp_bsln{1}); FCz_erp_bslnamp{2}.*cos(FCz_erp_bsln{2}); FCz_erp_bslnamp{3}.*cos(FCz_erp_bsln{3}); FCz_erp_bslnamp{4}.*cos(FCz_erp_bsln{4})];
        % erp baseline phase - sin(angle)
        continuous_variable(:, 4) = [FCz_erp_bslnamp{1}.*sin(FCz_erp_bsln{1}); FCz_erp_bslnamp{2}.*sin(FCz_erp_bsln{2}); FCz_erp_bslnamp{3}.*sin(FCz_erp_bsln{3}); FCz_erp_bslnamp{4}.*sin(FCz_erp_bsln{4})];
        % baseline angle
        continuous_variable(:, 5) = [FCz_erp_bsln{1}; FCz_erp_bsln{2}; FCz_erp_bsln{3}; FCz_erp_bsln{4}];
        % baseline amplitude
        continuous_variable(:, 6) = [FCz_erp_bsln_amplitude{1}'; FCz_erp_bsln_amplitude{2}'; FCz_erp_bsln_amplitude{3}'; FCz_erp_bsln_amplitude{4}'];
        
        % figure - plot effect of phase on amplitude
%         if  p == 37
%             figure; plot([FCz_erp_bsln{1}; FCz_erp_bsln{2}; FCz_erp_bsln{3}; FCz_erp_bsln{4}], ...
%                  mean(correct_go_trials_data(:, 241:240*1.5), 2), 'o', 'color', [.5 .5 .5])
%             ax = gca; c = ax.Color;
%             ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
%             xticks([-pi -pi/2 0 pi/2 pi])
%             xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
%             axis([-4 4 -inf inf]);
%             xlabel('Phase of slow fluctuations', 'FontSize', 32, 'FontWeight','normal')
%             ylabel('CNV amplitude \muV', 'FontSize', 32, 'FontWeight','normal')
%     %         title(['AB', num2str(p)], 'FontSize', 32, 'FontWeight','normal')
            
            % include only go/no-go task
%             figure; plot([FCz_erp_bsln{3}; FCz_erp_bsln{4}], ...
%              mean(correct_go_trials_data(length([FCz_erp_bsln{1}; FCz_erp_bsln{2}])+1:end, 241:240*1.5), 2), 'o', 'color', [0 0 0], 'MarkerSize', 10, 'LineWidth', 1.5)
%             ax = gca; c = ax.Color;
%             ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
%             xticks([-pi -pi/2 0 pi/2 pi])
%             xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
%             axis([-4 4 -inf inf]);
%             xlabel('Phase of slow fluctuations', 'FontSize', 32, 'FontWeight','normal')
%             ylabel('CNV amplitude \muV', 'FontSize', 32, 'FontWeight','normal')
%             title(['AB', num2str(p)], 'FontSize', 32, 'FontWeight','normal')
%         end
        
        
        % separate EEG baseline amplitude according to angle of baseline fluctuations
        data_sorted = sortrows([continuous_variable(:, 6), continuous_variable(:, 5)], 2);
        Bsln_amplitude_sorted{grp}(count, 1, :) = mean(data_sorted(1:floor(size(data_sorted, 1)/5), 1:end-1), 1);
        Bsln_amplitude_sorted{grp}(count, 2, :) = mean(data_sorted(floor(size(data_sorted, 1)/5)+1:2*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        Bsln_amplitude_sorted{grp}(count, 3, :) = mean(data_sorted(2*floor(size(data_sorted, 1)/5)+1:3*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        Bsln_amplitude_sorted{grp}(count, 4, :) = mean(data_sorted(3*floor(size(data_sorted, 1)/5)+1:4*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        Bsln_amplitude_sorted{grp}(count, 5, :) = mean(data_sorted(4*floor(size(data_sorted, 1)/5)+1:end, 1:end-1), 1);

         % separate ERP according to angle of baseline fluctuations
        data_sorted = sortrows([mean(correct_go_trials_data(:, 601:850), 2), continuous_variable(:, 5)], 2);
        CNV_angle_sorted{grp}(count, 1, :) = mean(data_sorted(1:floor(size(data_sorted, 1)/5), 1:end-1), 1);
        CNV_angle_sorted{grp}(count, 2, :) = mean(data_sorted(floor(size(data_sorted, 1)/5)+1:2*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        CNV_angle_sorted{grp}(count, 3, :) = mean(data_sorted(2*floor(size(data_sorted, 1)/5)+1:3*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        CNV_angle_sorted{grp}(count, 4, :) = mean(data_sorted(3*floor(size(data_sorted, 1)/5)+1:4*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        CNV_angle_sorted{grp}(count, 5, :) = mean(data_sorted(4*floor(size(data_sorted, 1)/5)+1:end, 1:end-1), 1);

        
        % create variable with avg data in 3 groups for each participant
        % effect of time on task, alpha, bsln cos, sin, angle, amplitude
        for var = 1:6
            data_sorted = sortrows([correct_go_trials_data, continuous_variable(:, var)], size(correct_go_trials_data, 2)+1);
            % CNV_tertiles = sorting variable, participant, tertile, time
            CNV_tertiles(grp, var, count, 1, :) = mean(data_sorted(1:floor(size(data_sorted, 1)/3), 1:end-1), 1);
            CNV_tertiles(grp, var, count, 2, :) = mean(data_sorted(floor(size(data_sorted, 1)/3)+1:2*floor(size(data_sorted, 1)/3), 1:end-1), 1);
            CNV_tertiles(grp, var, count, 3, :) = mean(data_sorted(2*floor(size(data_sorted, 1)/3)+1:end, 1:end-1), 1);
        end
        
        % create variable with run data for each participant
        CNV_run(grp, count, 1, :) = mean(correct_go_trials_data(categorical_variable(:,2) == 0, :), 1);
        CNV_run(grp, count, 2, :) = mean(correct_go_trials_data(categorical_variable(:,2) == 1, :), 1);  
    end
end

%% plot  CNV_angle_sorted

xaxis = -pi+2*pi/10*[1 3 5 7 9];

figure;
errorbar(xaxis, squeeze(mean(CNV_angle_sorted{1}, 1)), squeeze(std(CNV_angle_sorted{1}, [], 1))/sqrt(size(CNV_angle_sorted{1}, 1)), 'k', ...
    'LineWidth', 2)
hold on
errorbar(xaxis, squeeze(mean(CNV_angle_sorted{2}, 1)), squeeze(std(CNV_angle_sorted{2}, [], 1))/sqrt(size(CNV_angle_sorted{2}, 1)), '--r', ...
    'LineWidth', 2)

ax = gca; c = ax.Color;
ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
axis([-4 4 -inf inf]);
xlabel('Phase of slow fluctuations', 'FontSize', 32, 'FontWeight','normal')
ylabel('CNV amplitude \muV', 'FontSize', 32, 'FontWeight','normal')


%% plot  EEG bsln amplitude according to angle calculated in hilbert

xaxis = -pi+2*pi/10*[1 3 5 7 9];

figure;
errorbar(xaxis, squeeze(mean(Bsln_amplitude_sorted{1}, 1)), squeeze(std(Bsln_amplitude_sorted{1}, [], 1))/sqrt(size(Bsln_amplitude_sorted{1}, 1)), 'k', ...
    'LineWidth', 2)
hold on
errorbar(xaxis, squeeze(mean(Bsln_amplitude_sorted{2}, 1)), squeeze(std(Bsln_amplitude_sorted{2}, [], 1))/sqrt(size(Bsln_amplitude_sorted{2}, 1)), '--r', ...
    'LineWidth', 2)

ax = gca; c = ax.Color;
ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
axis([-4 4 -inf inf]);
xlabel('Phase of slow fluctuations', 'FontSize', 32, 'FontWeight','normal')
ylabel('Bsln EEG amp \muV', 'FontSize', 32, 'FontWeight','normal')


%% draw graphs
variable = {'Time-on-task','Alpha power' 'Cos(\theta)' 'Sin(\theta)', 'Phase at baseline', 'Amplitude at baseline'};
group_name = {'young', 'older'};

for var = 1:4
    for grp = 1:2
        clear data
        data = squeeze(CNV_tertiles(grp, var, :, :, 101:end)); % participant x tertile x time

        media1=squeeze(mean(data(:,1,:), 1));
        erropadrao1=squeeze(std(data(:,1,:), 0, 1)/sqrt(size(data(:,1,:), 1)));
        media2=squeeze(mean(data(:,2,:), 1));
        erropadrao2=squeeze(std(data(:,2,:), 0, 1)/sqrt(size(data(:,2,:), 1)));
        media3=squeeze(mean(data(:,3,:), 1));
        erropadrao3=squeeze(std(data(:,3,:), 0, 1)/sqrt(size(data(:,3,:), 1)));

        colormap cool;
        cmap = colormap;

        xaxis=(.002:1/500:6);
        figure;
        plot( xaxis', media1', 'color', cmap(1,:), 'linewidth', 2);
        hold on
        jbfill(xaxis,(media1+erropadrao1)',(media1-erropadrao1)', cmap(1,:), cmap(1,:), 1, 0.25)
        hold on
        plot( xaxis, media2, 'color', cmap(128,:), 'linewidth', 2);
        hold on
        jbfill(xaxis,(media2+erropadrao2)',(media2-erropadrao2)', cmap(128,:), cmap(128,:), 1, 0.25)
        hold on
        plot( xaxis', media3', 'color', cmap(256,:), 'linewidth', 2);
        hold on
        jbfill(xaxis,(media3+erropadrao2)',(media3-erropadrao2)', cmap(256,:), cmap(256,:), 1, 0.25)
        hold on
        plot(xaxis, zeros(1, size(data, 3)), ':k')
        hold off
        axis([-inf 3 -20 10]);
        ax = gca; c = ax.Color;
        ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        title([variable{var},' - ', group_name{grp}], 'FontSize', 32, 'FontWeight','normal')
    end
end

%% plot effect of run

for grp = 1:2
    clear data
    data = squeeze(CNV_run(grp,:, :, 101:end)); % participant x run x time

    media1=squeeze(mean(data(:,1,:), 1));
    erropadrao1=squeeze(std(data(:,1,:), 0, 1)/sqrt(size(data(:,1,:), 1)));
    media2=squeeze(mean(data(:,2,:), 1));
    erropadrao2=squeeze(std(data(:,2,:), 0, 1)/sqrt(size(data(:,2,:), 1)));

    colormap cool;
    cmap = colormap;

    xaxis=(.002:1/500:6);
    figure;
    plot( xaxis', media1', 'color', cmap(1,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(media1+erropadrao1)',(media1-erropadrao1)', cmap(1,:), cmap(1,:), 1, 0.25)
    hold on
    plot( xaxis, media2, 'color', cmap(256,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(media2+erropadrao2)',(media2-erropadrao2)', cmap(256,:), cmap(256,:), 1, 0.25)
    hold on
    plot(xaxis, zeros(1, size(data, 3)), ':k')
    hold off
    axis([-inf 3 -20 10]);
    ax = gca; c = ax.Color;
    ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    title(['Run - ', group_name{grp}], 'FontSize', 32, 'FontWeight','normal')
end