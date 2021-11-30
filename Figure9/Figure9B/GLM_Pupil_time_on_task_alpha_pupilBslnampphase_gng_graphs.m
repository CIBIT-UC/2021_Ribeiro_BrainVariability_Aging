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

betas_young = []; betas_older = [];
for grp = 1:2
    count = 0;
    for p = group{grp}
    
        count = count + 1;
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        % eeg directory
        eeg_dir = strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\EEG\');
        % pupil rejected epochs - these were manually/visually rejected
        pupil_directory=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\');
        load([pupil_directory, 'RejectEpochsCorrectGoCueOnly']); % data from G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\pupil_analysis_3b_RejectedEpochs_Go_CueOnly_SaveRunData.m
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
            EEG = pop_epoch( EEG, {  '1'  }, [-.2  6]); % epoch locked to cue onset
            % create variable with trial-by-trial pupil size at baseline as
            % percent signal change
%             pupil_baseline{t} = squeeze(mean(EEG.data(1, 1:48, setdiff(1:60, RejectEpochsCorrectGoCueOnly{t})), 2)); 
            EEG = pop_rmbase( EEG, [-200    0]); % subtract baseline pupil
            % delete epochs with artefacts as well as cue-only trials and error trials and trials after error trials
            % average across time from cue onset up to 1500 ms after the
            % cue (period between cue and target)
            pupil_data{t} = squeeze(EEG.data(1, 49:end, setdiff(1:60, RejectEpochsCorrectGoCueOnly{t}))); % channel x time x trials
            
%             % create variable with trial number for each run - time on task
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
        
        
       % calculate amplitude of pupil at cue-onset after filtering
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
                EEG = pop_epoch( EEG, {  '1'  }, [-.2  1.5]); % epoch locked to cue onset
                pupil_bsln_amplitude{t} = squeeze(EEG.data(1, 48, setdiff(1:60, RejectEpochsCorrectGoCueOnly{t}))); % amplitude at cue-onset for each trial
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
   
               
       %% match epochs of pupil, alpha, and erp
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
                pupil_data{t}(:, index2delete) = [];
                pupil_bsln_amplitude{t}(index2delete) = []; 
                time_on_task{t}(index2delete) = [];
            end
            % alpha
            index2delete = []; epochs2delete = setdiff(epochs_prestim_spectra{t}(:, 2), included_epochs{t});
            if ~isempty(epochs2delete)
                for epoch = 1:length(epochs2delete)
                    index2delete = [index2delete; find(epochs_prestim_spectra{t}(:, 2) == epochs2delete(epoch))];
                end
                alpha_power{t}(index2delete ) = []; 
%                 beta_power{t}(index2delete ) = []; 
%                 exponent_prestim{t}(index2delete ) = [];
            end
       end
       
            %% create variables for model
     correct_go_trials_data = [];  categorical_variable = []; 
     % pupil response: time x trials
        correct_go_trials_data=[pupil_data{1}'; pupil_data{2}'; pupil_data{3}'; pupil_data{4}']; % trials x time
        % categorical variable with task and run - one in each column
        categorical_variable(:,1)=[zeros(size(pupil_data{1}, 2)+size(pupil_data{2},2), 1); ones(size(pupil_data{3},2)+size(pupil_data{4},2), 1)];
        categorical_variable(:,2)=[zeros(size(pupil_data{1},2), 1); ones(size(pupil_data{2},2), 1);...
         zeros(size(pupil_data{3},2), 1); ones(size(pupil_data{4},2), 1)];

    % create continuous variable with time-on-task; baseline pupil ;alpha baseline - total 3 continuous variables
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
        continuous_variable(:, 4) = [pupil_bsln_amp{1}.*sin(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*sin(pupil_bsln_phase{2}); pupil_bsln_amp{3}.*sin(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*sin(pupil_bsln_phase{4})];
        % pupil baseline phase angle
        continuous_variable(:, 5) = [pupil_bsln_phase{1}; pupil_bsln_phase{2}; pupil_bsln_phase{3}; pupil_bsln_phase{4}];
        % pupil baseline amplitude
        continuous_variable(:, 6) = [pupil_bsln_amplitude{1}; pupil_bsln_amplitude{2}; pupil_bsln_amplitude{3}; pupil_bsln_amplitude{4}];
       
        
        if ismember(p, [34  44 45 56 62 72 74 76 78 81 19 22 73])
            % figure - plot effect of phase on amplitude
            figure; plot([pupil_bsln_phase{1}; pupil_bsln_phase{2}; pupil_bsln_phase{3}; pupil_bsln_phase{4}], ...
                 mean(correct_go_trials_data(:, 241:240*1.5), 2), 'o', 'color', [0 0 0], 'MarkerSize', 10, 'LineWidth', 1.5)
            ax = gca; c = ax.Color;
            ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
            xticks([-pi -pi/2 0 pi/2 pi])
            xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
            axis([-4 4 -inf inf]);
            xlabel('Phase of slow fluctuations', 'FontSize', 32, 'FontWeight','normal')
            ylabel('PD amplitude (%)', 'FontSize', 32, 'FontWeight','normal')
            title(['AB', num2str(p)], 'FontSize', 32, 'FontWeight','normal')
        end
        
        % separate pupil baseline amplitude according to angle of baseline fluctuations
        data_sorted = sortrows([continuous_variable(:, 6), continuous_variable(:, 5)], 2);
        Bsln_amplitude_sorted{grp}(count, 1, :) = mean(data_sorted(1:floor(size(data_sorted, 1)/5), 1:end-1), 1);
        Bsln_amplitude_sorted{grp}(count, 2, :) = mean(data_sorted(floor(size(data_sorted, 1)/5)+1:2*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        Bsln_amplitude_sorted{grp}(count, 3, :) = mean(data_sorted(2*floor(size(data_sorted, 1)/5)+1:3*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        Bsln_amplitude_sorted{grp}(count, 4, :) = mean(data_sorted(3*floor(size(data_sorted, 1)/5)+1:4*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        Bsln_amplitude_sorted{grp}(count, 5, :) = mean(data_sorted(4*floor(size(data_sorted, 1)/5)+1:end, 1:end-1), 1);

        
        % separate TPR according to angle of baseline fluctuations
        data_sorted = sortrows([mean(correct_go_trials_data(:, 241:240*1.5), 2), continuous_variable(:, 5)], 2);
        % pupil_tertiles = sorting variable, participant, tertile, time
        pupil_angle_sorted(grp, count, 1, :) = mean(data_sorted(1:floor(size(data_sorted, 1)/5), 1:end-1), 1);
        pupil_angle_sorted(grp, count, 2, :) = mean(data_sorted(floor(size(data_sorted, 1)/5)+1:2*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        pupil_angle_sorted(grp, count, 3, :) = mean(data_sorted(2*floor(size(data_sorted, 1)/5)+1:3*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        pupil_angle_sorted(grp, count, 4, :) = mean(data_sorted(3*floor(size(data_sorted, 1)/5)+1:4*floor(size(data_sorted, 1)/5), 1:end-1), 1);
        pupil_angle_sorted(grp, count, 5, :) = mean(data_sorted(4*floor(size(data_sorted, 1)/5)+1:end, 1:end-1), 1);
        
        % create variable with avg data in tertiles for each participant
        % effect of time on task, alpha, pupil baseline, eeg baseline
        for var = 1:4
            data_sorted = sortrows([correct_go_trials_data, continuous_variable(:, var)], size(correct_go_trials_data, 2)+1);
            % pupil_tertiles = sorting variable, participant, tertile, time
            pupil_tertiles(grp, var, count, 1, :) = mean(data_sorted(1:floor(size(data_sorted, 1)/3), 1:end-1), 1);
            pupil_tertiles(grp, var, count, 2, :) = mean(data_sorted(floor(size(data_sorted, 1)/3)+1:2*floor(size(data_sorted, 1)/3), 1:end-1), 1);
            pupil_tertiles(grp, var, count, 3, :) = mean(data_sorted(2*floor(size(data_sorted, 1)/3)+1:end, 1:end-1), 1);
        end
        
        % create variable with run data for each participant
        pupil_run(grp, count, 1, :) = mean(correct_go_trials_data(categorical_variable(:,2) == 0, :), 1);
        pupil_run(grp, count, 2, :) = mean(correct_go_trials_data(categorical_variable(:,2) == 1, :), 1);  
    end
end

%% plot  pupil_angle_sorted

xaxis = -pi+2*pi/10*[1 3 5 7 9];

figure;
errorbar(xaxis, squeeze(mean(pupil_angle_sorted(1, 1:35, :), 2)), squeeze(std(pupil_angle_sorted(1, 1:35, :), [], 2))/sqrt(35), 'k', ...
    'LineWidth', 2)
hold on
errorbar(xaxis, squeeze(mean(pupil_angle_sorted(2, :, :), 2)), squeeze(std(pupil_angle_sorted(2, :, :), [], 2))/sqrt(37), '--r', ...
    'LineWidth', 2)

ax = gca; c = ax.Color;
ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
axis([-4 4 -inf inf]);
xlabel('Phase of slow fluctuations', 'FontSize', 32, 'FontWeight','normal')
ylabel('PD amplitude (%)', 'FontSize', 32, 'FontWeight','normal')

%% plot  pupil bsln amplitude accroding to angle calculated in hilbert

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
ylabel('Bsln pupil amp (%)', 'FontSize', 32, 'FontWeight','normal')


%% draw graphs
variable = {'Time-on-task', 'Alpha power', 'Cos(\theta)', 'Sin(\theta)'};
group_name = {'young', 'older'};

for var = 1:4
    for grp = 1:2
        clear data
        data = squeeze(pupil_tertiles(grp, var, :, :, :)); % participant x tertile x time

        media1=squeeze(mean(data(:,1,:), 1));
        erropadrao1=squeeze(std(data(:,1,:), 0, 1)/sqrt(size(data(:,1,:), 1)));
        media2=squeeze(mean(data(:,2,:), 1));
        erropadrao2=squeeze(std(data(:,2,:), 0, 1)/sqrt(size(data(:,2,:), 1)));
        media3=squeeze(mean(data(:,3,:), 1));
        erropadrao3=squeeze(std(data(:,3,:), 0, 1)/sqrt(size(data(:,3,:), 1)));

        colormap cool;
        cmap = colormap;

        xaxis=0.001:1/240:6;
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
        axis([-inf inf -.06 .12]);
%         if task == 1
%             axis([-0.2 inf -.5 10]);
%         else
%             axis([-0.2 inf -.5 12]);
%         end
        ax = gca; c = ax.Color;
        ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Pupil diameter (%)', 'FontSize', 32, 'FontWeight','normal')
        title([variable{var},' - ', group_name{grp}], 'FontSize', 32, 'FontWeight','normal')
    end
end

%% plot effect of run

for grp = 1:2
    clear data
    data = squeeze(pupil_run(grp,:, :, :)); % participant x run x time

    media1=squeeze(mean(data(:,1,:), 1));
    erropadrao1=squeeze(std(data(:,1,:), 0, 1)/sqrt(size(data(:,1,:), 1)));
    media2=squeeze(mean(data(:,2,:), 1));
    erropadrao2=squeeze(std(data(:,2,:), 0, 1)/sqrt(size(data(:,2,:), 1)));

    colormap cool;
    cmap = colormap;

    xaxis=0.001:1/240:6;
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
    axis([-inf inf -.06 .12]);

    ax = gca; c = ax.Color;
    ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Pupil diameter (%)', 'FontSize', 32, 'FontWeight','normal')
    title(['Run - ', group_name{grp}], 'FontSize', 32, 'FontWeight','normal')
end
