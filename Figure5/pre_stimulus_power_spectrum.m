% calculation of pre-stimulus data spectra 1 sec before cue onset
% for calculation of pre-stimulus alpha amplitude - Maria Ribeiro 21May2020
clear; close all;

younger=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];
participants = {younger, older}; task={'W1', 'D1', 'D2', 'G1', 'G2'};
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for group = 2%1:2
    count = 0;
    for p = participants{group}
        count = count+1;
        participant=strcat('AB', num2str(p));
        dir_data=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\EEG\');
        for t = 1:5
            % clear eeglab
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            %load file - from eeglab_analysis_2b_synchronize_pupil
            % _RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica_ArtICsRem_Filt0_1_35Hz_Pupil.set
            % manual/visual removal of periods with artifacts
            % save as ..._Filt0_1_35Hz_Pupil_ManualArtRej.set
            if p == 77 && ismember(t, [2, 3])
                filename=strcat(participant, '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej_HEO.set');
            elseif ismember(p, [14, 15]) && t == 1
                filename=strcat(participant, '_W2_Filt0_1_35Hz_Pupil_ManualArtRej.set');
            else
                filename=strcat(participant, '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej.set');
            end
            EEG = pop_loadset(filename, dir_data);
            
            % interpolate bad channels removed earlier
            if length(EEG.chanlocs)<67 && p ~= 6 && p ~=38 % these participants do not have eyetracking data but has all channels ok
                [EEG] = interpol( EEG);
            end
            
            % epoch data - 3.5 sec before cue onset
            EEG = pop_epoch( EEG, {  '1'  }, [-3.5  0]);
            
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
            
            epochs_prestim = zeros(size(EEG.epoch, 2), 2); % epoch number x epoch original
            epochs_prestim(:, 1) = 1:size(EEG.epoch, 2);

            for epoch = 1:size(EEG.epoch, 2)
                for i = 1:size(EEG.event, 2)
                    if EEG.event(i).epoch == epoch && strcmp(EEG.event(i).type, '1')
                        epochs_prestim(epoch, 2) = trial_numbers((trial_numbers(:, 2) == EEG.event(i).urevent), 1);
                    end
                end
            end
            
            for channel = 1:59
%               EEG.data = channel, time, trials
                data4spectra = squeeze(EEG.data(channel, :, :));
                [Pxx, F] = pwelch(data4spectra, [], [], [], 500);% When X is a matrix, the PSD is
    %   computed independently for each column and stored in the corresponding
    %   column of Pxx. Note that the units of the PSD estimate are in squared magnitude units of the time series data per unit frequency.
    %    if the input data is in volts, the PSD estimate is in units of
    %    squared volts per unit frequency - https://www.mathworks.com/help/signal/ref/pwelch.html#btuf68p-3
                if group == 1
                    PowerSpectralDensity_Young{count, t, channel} = Pxx(1:50, :); % participant X task X electrode
                    epochs_prestim_young{count, t} = epochs_prestim;
                else
                    PowerSpectralDensity_Older{count, t, channel} = Pxx(1:50, :);
                    epochs_prestim_older{count, t} = epochs_prestim;
                end
%                 [spectra(channel),freqs,speccomp,contrib,specstd] =
%                 spectopo(data4spectra); % eeglab function
            end
        end
    end
end

Frequencies = F(1:50);
% save variables
save('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\PowerSpectralDensity_Young.mat', 'PowerSpectralDensity_Young');
save('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\PowerSpectralDensity_Older.mat', 'PowerSpectralDensity_Older');
save('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\Frequencies.mat', 'Frequencies');
save('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\epochs_prestim_older.mat', 'epochs_prestim_older');
save('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\epochs_prestim_young.mat', 'epochs_prestim_young');

%% plot data
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs.mat');

% plot power spectrum at FCz electrode
% find channel index
for c = 1:length(chanlocs)
    if strcmp(chanlocs(c).labels, 'FCz')
        chan = c;
    end
end

% plot spectrum simple RT task
for p = 1:size(PowerSpectralDensity_Young, 1)
    data_young_passive(p, :) = mean(PowerSpectralDensity_Young{p, 1, chan}, 2);
    data_young_simple_RT(p, :) = mean([PowerSpectralDensity_Young{p, 2, chan} PowerSpectralDensity_Young{p, 3, chan}], 2);
    data_young_gng(p, :) = mean([PowerSpectralDensity_Young{p, 4, chan} PowerSpectralDensity_Young{p, 5, chan}], 2);
end

for p = 1:size(PowerSpectralDensity_Older, 1)
    data_older_passive(p, :) = mean(PowerSpectralDensity_Older{p, 1, chan}, 2);
    data_older_simple_RT(p, :) = mean([PowerSpectralDensity_Older{p, 2, chan} PowerSpectralDensity_Older{p, 3, chan}], 2);
    data_older_gng(p, :) = mean([PowerSpectralDensity_Older{p, 4, chan} PowerSpectralDensity_Older{p, 5, chan}], 2);
end

% plot_spectrum_mean_std_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)
plot_spectrum_mean_std_young_older(data_young_passive, data_older_passive, Frequencies, 35, 'Passive', 'Power (dB)')
plot_spectrum_mean_std_young_older(data_young_simple_RT, data_older_simple_RT, Frequencies, 35, 'Simple RT', 'Power (dB)')
plot_spectrum_mean_std_young_older(data_young_gng, data_older_gng, Frequencies, 35, 'Go/no-go', 'Power (dB)')

%% plot graph function - group comparison
function plot_spectrum_mean_std_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)

    mean_data_young=mean(data_young(:, frequencies <= max_freq2plot), 1);
    se_data_young=std(data_young(:, frequencies <= max_freq2plot), 0, 1)/sqrt(size(data_young, 1));

    mean_data_older=mean(data_older(:, frequencies <= max_freq2plot), 1);
    se_data_older=std(data_older(:, frequencies <= max_freq2plot), 0, 1)/sqrt(size(data_older, 1));


    % color gng young = [0 .5 0]
    % color simple RT young =   [1 .5 0]
    % color gng older = [0 0 .5]
    % color simple RT older =  [.75 0.25 0]

    x_axis = frequencies(frequencies <= max_freq2plot); % frequencies
    figure;
    % % plot a line at zero
    % plot([0 0],[0 11], '--', 'color', [0 0 0]);
    % hold on
    % detection data
    loglog(x_axis, mean_data_young,'color', 'k', 'LineWidth',1.5)
    hold on
    jbfill(x_axis',mean_data_young+se_data_young, mean_data_young-se_data_young, 'k','k', 0.1)
    % GNG data
    hold on
    loglog(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
    hold on
    jbfill(x_axis', mean_data_older+se_data_older, mean_data_older-se_data_older,'r','r', 0.1)
    hold off
    ax = gca;
    % c = ax.Color;
    % legend('Detection', 'GNG')
    ax.FontSize = 20;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    axis([1 30 -inf inf])
    xlabel('Frequency (Hz)', 'FontSize', 28, 'FontWeight','normal')
    ylabel(y_axis_label, 'FontSize', 28, 'FontWeight','normal')
    title(graph_title, 'FontSize', 28, 'FontWeight','normal')

end


%% function to interpolate missing channels - Marco Simões   
% interpolates and reorders channels so they are in right order!
    function [EEG] = interpol( EEG, chanlocs )
        % interpolation
        if nargin < 2
            load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs.mat');
        end

        chans_eeg = [];
        for i=1:length(EEG.chanlocs)
            chans_eeg = [ chans_eeg {EEG.chanlocs(i).labels} ];
        end

        idxs = [];
        for i=1:length(chanlocs)
            index = find(ismember(chans_eeg, chanlocs(i).labels) == 1, 1);
            if isempty(index)
                idxs = [idxs i];
            end
        end

        EEG = pop_interp(EEG, chanlocs(idxs), 'spherical');

        % reorder
        chans_eeg = [];
        for c=1:length(EEG.chanlocs)
            chans_eeg = [ chans_eeg {EEG.chanlocs(c).labels} ];
        end

        idxs = [];
        for c=1:length(chanlocs)
            index = find(ismember(chans_eeg, chanlocs(c).labels) == 1, 1);
            idxs = [idxs index];
        end

        %if length(idxs) == 58

           EEG.data = EEG.data(idxs,:,:);
           EEG.chanlocs = EEG.chanlocs(idxs);

           indcomps = [];
           for compidx = 1:length(EEG.icachansind)
               indcomps = [indcomps find(EEG.icachansind(compidx) == idxs)];
           end
           EEG.icachansind = indcomps;

   % end
end
            