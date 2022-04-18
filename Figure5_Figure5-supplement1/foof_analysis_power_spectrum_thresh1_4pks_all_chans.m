% fooof analysis at chosen electrode - pre-stimulus epoch 3.5 sec before cue
% onset
clear; close all

% load data - data were calculated here:
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\pre_stimulus_power_spectrum.m
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra';
load PowerSpectralDensity_Young; load PowerSpectralDensity_Older; load Frequencies;

load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs.mat');

% loop through electrodes
for chan = 1:59
    % avg power spectrum across trials
    % group: 1 = young; 2 = older; task: 1 = passive; 2 = simple RT; 3 = gng
    for p = 1:size(PowerSpectralDensity_Young, 1)
        avg_power_spectrum_young(1, p, :) = mean(PowerSpectralDensity_Young{p, 1, chan}, 2);
        avg_power_spectrum_young(2, p, :) = mean([PowerSpectralDensity_Young{p, 2, chan} PowerSpectralDensity_Young{p, 3, chan}], 2);
        avg_power_spectrum_young(3, p, :) = mean([PowerSpectralDensity_Young{p, 4, chan} PowerSpectralDensity_Young{p, 5, chan}], 2);
    end

    for p = 1:size(PowerSpectralDensity_Older, 1)
        avg_power_spectrum_older(1, p, :) = mean(PowerSpectralDensity_Older{p, 1, chan}, 2);
        avg_power_spectrum_older(2, p, :) =  mean([PowerSpectralDensity_Older{p, 2, chan} PowerSpectralDensity_Older{p, 3, chan}], 2);
        avg_power_spectrum_older(3, p, :) = mean([PowerSpectralDensity_Older{p, 4, chan} PowerSpectralDensity_Older{p, 5, chan}], 2);
    end

    %% fooof analyses
    % best settings to avoid overfitting:
    % https://fooof-tools.github.io/fooof/auto_tutorials/plot_07-TroubleShooting.html#sphx-glr-auto-tutorials-plot-07-troubleshooting-py
    % 
    %   settings        = struct, can optionally include:
    %       settings.peak_width_limits
    %       settings.max_n_peaks
    %       settings.min_peak_height
    %       settings.peak_threshold
    %       settings.aperiodic_mode
    %       settings.verbose
    % settings.peak_width_limits = [1 15]; settings.min_peak_height = 0.1;
    % aperiodic_mode = ‘fixed’, peak_width_limits = [2, inf], peak_threshold = 1, and default
    % settings otherwise). from 2020 Tran et al -JCognNeuroscii Linked sources of neural noise contribute to age- related cognitive decline
    settings.peak_width_limits = [1 8]; settings.peak_threshold = 1; settings.aperiodic_mode = 'fixed'; settings.min_peak_height = 0.1;
    settings.max_n_peaks = 4;
    % young
    for task = 1:3
        for p = 1:size(avg_power_spectrum_young, 2)
            count_theta = 0; count_alpha = 0; count_beta = 0;
            % fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
            fooof_results = fooof(Frequencies', squeeze(avg_power_spectrum_young(task, p, :))', [1 35], settings, 1);
    %         fooof_plot(fooof_results)
            exponent_young(task, p, chan) = fooof_results.aperiodic_params(2);
            offset_young(task, p, chan) = fooof_results.aperiodic_params(1);
            oscillations_power_young(task, p, chan, :) = [0 0 0]; % theta alpha beta - if oscillations not detected power = 0;
            oscillations_frequency_young(task, p, chan, :) = [NaN NaN NaN];
            oscillations_peak_width_young(task, p, chan, :) = [NaN NaN NaN];
            for peak = 1:size(fooof_results.peak_params, 1)
                oscillations_young{p}(task, chan, peak, :) = fooof_results.peak_params(peak, :);
                if fooof_results.peak_params(peak, 1)>=4 && fooof_results.peak_params(peak, 1)<=7  % theta - previous 8
                    count_theta = count_theta + 1;
                    % if there are two peaks within the same frequency band,use
                    % peak with highest power
                    if fooof_results.peak_params(peak, 2) > oscillations_power_young(task, p, 1)
                        oscillations_power_young(task, p, chan, 1) = fooof_results.peak_params(peak, 2); % power
                        oscillations_frequency_young(task, p, chan, 1) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width_young(task, p, chan, 1) = fooof_results.peak_params(peak, 3); % band width
                    end
                elseif fooof_results.peak_params(peak, 1)>7 && fooof_results.peak_params(peak, 1)< 14 % alpha previous 8 - 13
                    count_alpha = count_alpha + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power_young(task, p, chan, 2)
                        oscillations_power_young(task, p, chan, 2) = fooof_results.peak_params(peak, 2);
                        oscillations_frequency_young(task, p, chan, 2) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width_young(task, p, chan, 2) = fooof_results.peak_params(peak, 3); % band width
                    end
                elseif fooof_results.peak_params(peak, 1)>=14 && fooof_results.peak_params(peak, 1)< 30 % beta
                    count_beta = count_beta + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power_young(task, p, chan, 3)
                        oscillations_power_young(task, p, chan, 3) = fooof_results.peak_params(peak, 2);
                        oscillations_frequency_young(task, p, chan, 3) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width_young(task, p, chan, 3) = fooof_results.peak_params(peak, 3); % band width
                    end
                end
            end
            r_squared_error_young(task, p, chan, :) = [fooof_results.r_squared fooof_results.error];
            number_peaks_young(task, p, chan, :) = [count_theta; count_alpha; count_beta];
        end
    end

    % older
    for task = 1:3
        for p = 1:size(avg_power_spectrum_older, 2)
            count_theta = 0; count_alpha = 0; count_beta = 0;
            % fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
            fooof_results = fooof(Frequencies', squeeze(avg_power_spectrum_older(task, p, :))', [1 35], settings, 1);
    %         fooof_plot(fooof_results)
            exponent_older(task, p, chan) = fooof_results.aperiodic_params(2);
            offset_older(task, p, chan) = fooof_results.aperiodic_params(1);
            oscillations_power_older(task, p, chan, :) = [0 0 0]; % theta alpha beta - if oscillations not detected power = 0;
            oscillations_frequency_older(task, p, chan, :) = [NaN NaN NaN];
            oscillations_peak_width_older(task, p, chan, :) = [NaN NaN NaN];
            for peak = 1:size(fooof_results.peak_params, 1)
                oscillations_older{p}(task, chan, peak, :) = fooof_results.peak_params(peak, :);
                if fooof_results.peak_params(peak, 1)>=4 && fooof_results.peak_params(peak, 1)<=7  % theta - previous 8
                    count_theta = count_theta + 1;
                    % if there are two peaks within the same frequency band,use
                    % peak with highest power
                    if fooof_results.peak_params(peak, 2) > oscillations_power_older(task, p, 1)
                        oscillations_power_older(task, p, chan, 1) = fooof_results.peak_params(peak, 2); % power
                        oscillations_frequency_older(task, p, chan, 1) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width_older(task, p, chan, 1) = fooof_results.peak_params(peak, 3); % band width
                    end
                elseif fooof_results.peak_params(peak, 1)>7 && fooof_results.peak_params(peak, 1)< 14 % alpha previous 8 - 13
                    count_alpha = count_alpha + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power_older(task, p, chan, 2)
                        oscillations_power_older(task, p, chan, 2) = fooof_results.peak_params(peak, 2);
                        oscillations_frequency_older(task, p, chan, 2) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width_older(task, p, chan, 2) = fooof_results.peak_params(peak, 3); % band width
                    end
                elseif fooof_results.peak_params(peak, 1)>=14 && fooof_results.peak_params(peak, 1)< 30 % beta
                    count_beta = count_beta + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power_older(task, p, chan, 3)
                        oscillations_power_older(task, p, chan, 3) = fooof_results.peak_params(peak, 2);
                        oscillations_frequency_older(task, p, chan, 3) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width_older(task, p, chan, 3) = fooof_results.peak_params(peak, 3); % band width
                    end
                end
            end
            r_squared_error_older(task, p, chan, :) = [fooof_results.r_squared fooof_results.error];
            number_peaks_older(task, p, chan, :) = [count_theta; count_alpha; count_beta];
        end
    end
end
%% save results
save_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
cd(save_dir)
save exponent_young exponent_young
save offset_young offset_young
save oscillations_power_young oscillations_power_young 
save oscillations_frequency_young oscillations_frequency_young
save oscillations_peak_width_young oscillations_peak_width_young 
save oscillations_young oscillations_young
save r_squared_error_young r_squared_error_young 
save number_peaks_young number_peaks_young

save exponent_older exponent_older
save offset_older offset_older
save oscillations_power_older oscillations_power_older 
save oscillations_frequency_older oscillations_frequency_older
save oscillations_peak_width_older oscillations_peak_width_older 
save oscillations_older oscillations_older
save r_squared_error_older r_squared_error_older 
save number_peaks_older number_peaks_older

%% exclude outliers where the fittign wasn't so good - for each channel find participants where R2 zscore > 2.5
% change data for NAN - data in paper including all data points - outliers
% opted at the end for not excluding these outliers - results qualitatively
% the same
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
load r_squared_error_older % (task, p, chan, r2 and error)
load r_squared_error_young
load exponent_young; load exponent_older; load offset_young; load offset_older; 
load oscillations_power_young; load oscillations_power_older

% find participants where average R2 across channels was not so good -
% exclude those participants from analyses with fooof parameters
% for task = 1:3
%     for chan = 1:59
%         R2_transf_young(task, chan, :) = squeeze(fisherz(sqrt(r_squared_error_young(task, :, chan, 1))));
%     end
%     outliers_young{task} = find(abs(zscore(squeeze(mean(R2_transf_young(task, :, :), 2)))) > 2.5); % outliers
% end
% 
%     incl_young{task, chan} = find(abs(zscore(R2_transf_young)) <= 2.5); % without outliers
%     %exponent_young(task, p, chan)
%     exponent_young(task, outliers_young{task, chan}, chan) = NaN;
%     offset_young(task, outliers_young{task, chan}, chan) = NaN;
%     oscillations_power_young(task, outliers_young{task, chan}, chan, :) = NaN;
% 
% 
%     R2_transf_older = squeeze(fisherz(sqrt(r_squared_error_older(task, :, chan, 1))));
%     outliers_older{task, chan} = find(abs(zscore(R2_transf_older)) > 2.5); % outliers
%     incl_older{task, chan} = find(abs(zscore(R2_transf_older)) <= 2.5);
%     exponent_older(task, outliers_older{task, chan}, chan) = NaN;
%     offset_older(task, outliers_older{task, chan}, chan) = NaN;
%     oscillations_power_older(task, outliers_older{task, chan}, chan, :) = NaN;

% end


for task = 1:3
    for chan = 1:59
        R2_transf_young(:, chan) = squeeze(fisherz(sqrt(r_squared_error_young(task, :, chan, 1))));
%         outliers_young{task, chan} = find(abs(zscore(R2_transf_young)) > 3); % outliers
%         incl_young{task, chan} = find(abs(zscore(R2_transf_young)) <= 2.5); % without outliers
%         %exponent_young(task, p, chan)
%         exponent_young(task, outliers_young{task, chan}, chan) = NaN;
%         offset_young(task, outliers_young{task, chan}, chan) = NaN;
%         oscillations_power_young(task, outliers_young{task, chan}, chan, :) = NaN;
        
        
        R2_transf_older(:, chan) = squeeze(fisherz(sqrt(r_squared_error_older(task, :, chan, 1))));
%         outliers_older{task, chan} = find(abs(zscore(R2_transf_older)) > 3); % outliers
%         incl_older{task, chan} = find(abs(zscore(R2_transf_older)) <= 2.5);
%         exponent_older(task, outliers_older{task, chan}, chan) = NaN;
%         offset_older(task, outliers_older{task, chan}, chan) = NaN;
%         oscillations_power_older(task, outliers_older{task, chan}, chan, :) = NaN;
%         
    end
    
    zscoreR2_young(task, :, :) = zscore(R2_transf_young, 0, 'all');
    zscoreR2_older(task, :, :) = zscore(R2_transf_older, 0, 'all');
    [index_part_older{task}, index_chan_older{task}] = find(abs(zscore(R2_transf_older)) > 3);
    outliers_older{task} = unique(index_part_older{task});
    [index_part_young{task}, index_chan_young{task}] = find(abs(zscore(R2_transf_young)) > 3);
    outliers_young{task} = unique(index_part_young{task});
    
end

%% display topography of various variables
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
load exponent_young; load exponent_older;
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat

title_txt = {'Passive' 'Simple RT' 'Go/no-go'};

for task = 3%1:3

%         exponent_young(task, p, chan)
    figure;
    topoplot(mean(exponent_young(task, :, :), 2,'omitnan'), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([0 1.5]); %c.Axis.FontSize = 16;
%     colorbar;
%     colorbar('Ticks',[0, .5, 1, 1.5], 'FontSize', 24, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
%     title(['Young - ', title_txt{task}], 'FontSize', 24, 'FontWeight','normal')
%     set(get(gca,'title'),'Position',[0,-.8, 0])
%     text(5, 0.4, title_txt{task, grp})
    
    figure;
    topoplot(mean(exponent_older(task, :, :), 2,'omitnan'), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([0 1.5]); %c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[0.0, .5, 1.0, 1.5], 'FontSize', 40, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
%     title(['Older - ', title_txt{task}], 'FontSize', 24, 'FontWeight','normal')
%     set(get(gca,'title'),'Position',[0,-.8, 0])
%     text(5, 0.4, title_txt{task, grp})

end


%% stats using t-test with significance level corrected for bonferoni - 59 comparisons for each electrode
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
load exponent_young; load exponent_older;
for task = 3%1:3
    
    % Cite As
    % David Groppe (2021). mult_comp_perm_t2(data1,data2,n_perm,tail,alpha_level,mu,t_stat,reports,seed_state) (https://www.mathworks.com/matlabcentral/fileexchange/54585-mult_comp_perm_t2-data1-data2-n_perm-tail-alpha_level-mu-t_stat-reports-seed_state), MATLAB Central File Exchange. Retrieved April 16, 2021.
     [pval_t2(task, :), t_orig_t2(task, :), crit_t_t2(task, :), est_alpha, seed_state] = mult_comp_perm_t2(squeeze(exponent_young(task, :, :)), squeeze(exponent_older(task, :, :)));

%     figure;
    sig_t_values = zeros(1, 59);
    sig_t_values(t_orig_t2(task, :) > crit_t_t2(task, 2)) = 1;
    sig_t_values(t_orig_t2(task, :) < crit_t_t2(task, 1)) = 1;
    sig_chan_number = find(sig_t_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end


    figure;
    topoplot(ones(59, 1)*.030, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
        'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
    
    topoplot(t_orig_t2(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([0 8]); %c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[0,2, 4, 6, 8], 'FontSize', 40, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
%     title(['Young - ', title_txt{task}], 'FontSize', 24, 'FontWeight','normal')
%     set(get(gca,'title'),'Position',[0,-.8, 0])
%     text(5, 0.4, title_txt{task, grp})
    
end


%% plot offset topographies

load offset_young; load offset_older; 
for task = 1:3

%         exponent_young(task, p, chan)
    figure;
    topoplot(mean(offset_young(task, :, :), 2,'omitnan'), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([0 1.5]); %c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[0, .5, 1, 1.5], 'FontSize', 24, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
    
    figure;
    topoplot(mean(offset_older(task, :, :), 2,'omitnan'), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([0 1.5]); %c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[0, .5, 1, 1.5], 'FontSize', 24, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php

end

%% stats permutation
load offset_young; load offset_older; 
for task = 3%1:3
     [pval_t2(task, :), t_orig_t2(task, :), crit_t_t2(task, :), est_alpha, seed_state] = mult_comp_perm_t2(squeeze(offset_young(task, :, :)), squeeze(offset_older(task, :, :)));

%     figure;
    sig_t_values = zeros(1, 59);
    sig_t_values(t_orig_t2(task, :) > crit_t_t2(task, 2)) = 1;
    sig_t_values(t_orig_t2(task, :) < crit_t_t2(task, 1)) = 1;
    sig_chan_number = find(sig_t_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end

    figure;
    topoplot(ones(59, 1)*.03, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
        'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
    
    topoplot(t_orig_t2(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([0 8]); %c.Axis.FontSize = 16;
%     colorbar;
%     colorbar('Ticks',[0,2, 4, 6, 8], 'FontSize', 48, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php

end

%% plot topography of alpha and beta amplitude

load oscillations_power_young; load oscillations_power_older
% oscillations_power_older(task, p, chan, freq band) 

for task = 3%2:3
    for freq_band = 2

        figure;
        topoplot(mean(oscillations_power_young(task, :, :, freq_band), 2,'omitnan'), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        if freq_band == 2
            caxis([0 1]); 
        else
            caxis([0 .5]); %c.Axis.FontSize = 16;
        end
        colorbar;
        if freq_band == 2
            colorbar('Ticks',[0, .5, 1], 'FontSize',  40, 'FontWeight','normal');
        else
            colorbar('Ticks',[0, .25, .5], 'FontSize',  40, 'FontWeight','normal');
        end
        colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php

        figure;
        topoplot(mean(oscillations_power_older(task, :, :, freq_band), 2,'omitnan'), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        if freq_band == 2
            caxis([0 1]); 
        else
            caxis([0 .5]); %c.Axis.FontSize = 16;
        end
        colorbar;
        if freq_band == 2
            colorbar('Ticks',[0, .5, 1], 'FontSize',  40, 'FontWeight','normal');
        else
            colorbar('Ticks',[0, .25, .5], 'FontSize',  40, 'FontWeight','normal');
        end
        colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php

    end

end

%% stats 
load oscillations_power_young; load oscillations_power_older % not excluding outliers
for freq_band = 3 % band 1 = theta; 2 = alpha; 3 = beta
    for task = 3%1:3        
         [pval_t2(task, :), t_orig_t2(task, :), crit_t_t2(task, :), est_alpha, seed_state] = mult_comp_perm_t2(squeeze(oscillations_power_young(task, :, :, freq_band)), squeeze(oscillations_power_older(task, :, :, freq_band)));

%         figure;
        sig_t_values = zeros(1, 59);
        sig_t_values(t_orig_t2(task, :) > crit_t_t2(task, 2)) = 1;
        sig_t_values(t_orig_t2(task, :) < crit_t_t2(task, 1)) = 1;
        sig_chan_number = find(sig_t_values == 1);

        for x = 1:length(sig_chan_number)
            sig_chans{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
        end

        figure;
        topoplot(ones(59, 1)*.03, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
            'plotdisk', 'on',  'hcolor'  , 'none') ; hold on

        topoplot(t_orig_t2(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        colorbar;
        if freq_band == 3
            caxis([-5 0]); %c.Axis.FontSize = 16;
            colorbar('Ticks',[-4, -2, 0], 'FontSize', 40, 'FontWeight','normal');
            colormap(flip(crameri('imola'), 1)); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
        else
            caxis([0 3]); %c.Axis.FontSize = 16;
            colorbar('Ticks',[0, 1, 2, 3], 'FontSize', 40, 'FontWeight','normal');
            colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
        end

    end

end

%% compare parameters across tasks
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
load exponent_young; load exponent_older; load offset_young; load offset_older; 
load oscillations_power_young; load oscillations_power_older
% [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(data,n_perm,tail,alpha_level,mu,reports,seed_state)
[pval_t1(1, :), t_orig_t1(1, :), crit_t_t1(1, :), est_alpha, seed_state] = mult_comp_perm_t1(squeeze(exponent_young(2, :, :))-squeeze(exponent_young(3, :, :)));
[pval_t1(2, :), t_orig_t1(2, :), crit_t_t1(2, :), est_alpha, seed_state] = mult_comp_perm_t1(squeeze(exponent_older(2, :, :))-squeeze(exponent_older(3, :, :)));


[pval_t1(1, :), t_orig_t1(1, :), crit_t_t1(1, :), est_alpha, seed_state] = mult_comp_perm_t1(squeeze(offset_young(2, :, :))-squeeze(offset_young(3, :, :)));
[pval_t1(2, :), t_orig_t1(2, :), crit_t_t1(2, :), est_alpha, seed_state] = mult_comp_perm_t1(squeeze(offset_older(2, :, :))-squeeze(offset_older(3, :, :)));

% alpha power
[pval_t1(1, :), t_orig_t1(1, :), crit_t_t1(1, :), est_alpha, seed_state] = mult_comp_perm_t1(squeeze(oscillations_power_young(2, :, :, 2))-squeeze(oscillations_power_young(3, :, :, 2)));
[pval_t1(2, :), t_orig_t1(2, :), crit_t_t1(2, :), est_alpha, seed_state] = mult_comp_perm_t1(squeeze(oscillations_power_older(2, :, :, 2))-squeeze(oscillations_power_older(3, :, :, 2)));

% beta power
[pval_t1(1, :), t_orig_t1(1, :), crit_t_t1(1, :), est_alpha, seed_state] = mult_comp_perm_t1(squeeze(oscillations_power_young(2, :, :, 3))-squeeze(oscillations_power_young(3, :, :, 3)));
[pval_t1(2, :), t_orig_t1(2, :), crit_t_t1(2, :), est_alpha, seed_state] = mult_comp_perm_t1(squeeze(oscillations_power_older(2, :, :, 3))-squeeze(oscillations_power_older(3, :, :, 3)));

%% active tasks - simple RT vs gng
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
% load r_squared_error_older % (task, p, chan, r2 and error)
% load r_squared_error_young
load exponent_young; load exponent_older; load offset_young; load offset_older; 
load oscillations_power_young; load oscillations_power_older
%  exponent_young(task, p, chan)
% load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat

% plot FCz channel = 16
chan = 16;
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
% plot_all_data_2tasks(exponent_young(2, :, chan), exponent_young(3, :, chan), exponent_older(2, :, chan), exponent_older(3, :, chan), 'Exponent')
% plot_all_data_2tasks(offset_young(2, :, chan), offset_young(3, :, chan), offset_older(2, :, chan), offset_older(3, :, chan), 'Offset')
% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(exponent_young(3, :, chan), exponent_older(3, :, chan), 'Exponent')
plot_all_data_onetask(offset_young(3, :, chan), offset_older(3, :, chan), 'Offset')
%
% change zeros in power to NaN for graphs and SPSS stats
% oscillations_power_older(task, p, chan, freq band) 
oscillations_power_young(oscillations_power_young == 0) = NaN;
oscillations_power_older(oscillations_power_older == 0) = NaN;

y_axis_text = {'Theta power' 'Alpha power', 'Beta power'};
% y_axis_text_freq = {'Theta frequency (Hz)' 'Alpha frequency(Hz)', 'Beta frequency (Hz)'};

for peak = 2:3
%     plot_all_data_2groups(squeeze(oscillations_power_young(1, :, chan, peak)), squeeze(oscillations_power_older(1, :, chan, peak)),y_axis_text{peak});
%     plot_all_data_2tasks(squeeze(oscillations_power_young(2, :, chan, peak )),squeeze(oscillations_power_young(3, :, chan, peak )),...
%         squeeze(oscillations_power_older(2, :, chan, peak )), squeeze(oscillations_power_older(3, :, chan, peak )),y_axis_text{peak})
    
% % plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(squeeze(oscillations_power_young(3, :, chan, peak ))', squeeze(oscillations_power_older(3, :, chan, peak ))',y_axis_text{peak})

end


%% plot_spectrum_mean_std_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra';
load PowerSpectralDensity_Young; load PowerSpectralDensity_Older; load Frequencies;

load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_text = {'Passive', 'Simple RT', 'Go/no-go'};

% spectrum for FCz (chan 16) and POZ (chan 51)
for chan = [16, 51]
    clear avg_power_spectrum_young avg_power_spectrum_older
    % avg power spectrum across trials
    % group: 1 = young; 2 = older; task: 1 = passive; 2 = simple RT; 3 = gng
    for p = 1:size(PowerSpectralDensity_Young, 1)
        avg_power_spectrum_young(1, p, :) = mean(PowerSpectralDensity_Young{p, 1, chan}, 2);
        avg_power_spectrum_young(2, p, :) = mean([PowerSpectralDensity_Young{p, 2, chan} PowerSpectralDensity_Young{p, 3, chan}], 2);
        avg_power_spectrum_young(3, p, :) = mean([PowerSpectralDensity_Young{p, 4, chan} PowerSpectralDensity_Young{p, 5, chan}], 2);
    end

    for p = 1:size(PowerSpectralDensity_Older, 1)
        avg_power_spectrum_older(1, p, :) = mean(PowerSpectralDensity_Older{p, 1, chan}, 2);
        avg_power_spectrum_older(2, p, :) =  mean([PowerSpectralDensity_Older{p, 2, chan} PowerSpectralDensity_Older{p, 3, chan}], 2);
        avg_power_spectrum_older(3, p, :) = mean([PowerSpectralDensity_Older{p, 4, chan} PowerSpectralDensity_Older{p, 5, chan}], 2);
    end


    for task = 3%1:3
%         plot_spectrum_mean_std_young_older(log10(squeeze(avg_power_spectrum_young(task, :, :))), log10(squeeze(avg_power_spectrum_older(task, :, :))), Frequencies, 35, [chanlocs_EEGChanOnly(chan).labels, ' - ' title_text{task}], 'log10(PSD)')
        plot_spectrum_mean_std_young_older(log10(squeeze(avg_power_spectrum_young(task, :, :))), log10(squeeze(avg_power_spectrum_older(task, :, :))), Frequencies, 35, [chanlocs_EEGChanOnly(chan).labels], 'log10(PSD)')
    end
end

%% R2_square
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
load r_squared_error_older % (task, p, chan, r2 and error)
load r_squared_error_young


load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat

for task = 1:3
    R2_transf_young = []; R2_transf_older = [];
    for chan = 1:59
        R2_transf_young(:, chan) = squeeze(fisherz(sqrt(r_squared_error_young(task, :, chan, 1))));
        R2_transf_older(:, chan) = squeeze(fisherz(sqrt(r_squared_error_older(task, :, chan, 1))));
    end
    
    [pval_t2(task, :), t_orig_t2(task, :), crit_t_t2(task, :), est_alpha, seed_state] = mult_comp_perm_t2(R2_transf_young, R2_transf_older);

    sig_t_values = zeros(1, 59);
    sig_t_values(t_orig_t2(task, :) > crit_t_t2(task, 2)) = 1;
    sig_t_values(t_orig_t2(task, :) < crit_t_t2(task, 1)) = 1;
    sig_chan_number = find(sig_t_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end

    figure;
    topoplot(ones(59, 1)*.03, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
        'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
    
    topoplot(t_orig_t2(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([0 6]); %c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[0,2, 4, 6], 'FontSize', 48, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php

end


%% plot R2 at FCz = chan 16

R2_transf_young = []; R2_transf_older = []; chan = 16;
for task = 1:3
    R2_transf_young(task, :) = squeeze(fisherz(sqrt(r_squared_error_young(task, :, chan, 1))));
    R2_transf_older(task, :) = squeeze(fisherz(sqrt(r_squared_error_older(task, :, chan, 1))));
end

plot_all_data_2groups(R2_transf_young(1, :), R2_transf_older(1, :),'Fisher r-z');

plot_all_data_2tasks(R2_transf_young(2, :), R2_transf_young(3, :),...
    R2_transf_older(2, :), R2_transf_older(3, :), 'FisherZ(R)')

% plot only go/no-go task
% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(R2_transf_young(3, :), R2_transf_older(3, :), {'{\itR} (z-transformed)'})


%% across participants average R2 (average across electrodes)
mean(mean(r_squared_error_young(2, :, :, 1), 3), 2)
mean(mean(r_squared_error_young(3, :, :, 1), 3), 2)

mean(mean(r_squared_error_older(2, :, :, 1), 3), 2)
mean(mean(r_squared_error_older(3, :, :, 1), 3), 2)


%% error
% plot_all_data_points(r_squared_error_young(1, :, 2), r_squared_error_young(2, :, 2), r_squared_error_young(3, :, 2),...
%     r_squared_error_older(1, :, 2), r_squared_error_older(2, :, 2), r_squared_error_older(3, :, 2), 'Error');

plot_all_data_2groups(r_squared_error_young(1, :, 2), r_squared_error_older(1, :, 2),'Error');

plot_all_data_2tasks( r_squared_error_young(2, :, 2), r_squared_error_young(3, :, 2), r_squared_error_older(2, :, 2), r_squared_error_older(3, :, 2),  'Error')

plot_all_data_2tasks( fisherz(sqrt(r_squared_error_young(2, :, 1))), fisherz(sqrt(r_squared_error_young(3, :, 1))), ...
    fisherz(sqrt(r_squared_error_older(2, :, 1))), fisherz(sqrt(r_squared_error_older(3, :, 1))),  'Fisher r to z')


%% scatter plots - R2 (transformed) vs exponent and offset
load exponent_young; load exponent_older; load offset_young; load offset_older; 
%young
task_name = {'Simple RT', 'Go/no-go'};
for task=1:2
    R2_transf = squeeze(fisherz(sqrt(r_squared_error_young(task+1, :, 1))));
    outliers = find(abs(zscore(R2_transf)) > 2.5);
    figure;
    plot(R2_transf(abs(zscore(R2_transf)) <= 2.5),exponent_young(task+1, abs(zscore(R2_transf)) <= 2.5), 'o', 'color', 'k',  'MarkerSize',10); lsline
    hold on; plot(R2_transf(outliers),exponent_young(task+1, outliers), 'o', 'color', 'k', 'MarkerFaceColor', 'r',  'MarkerSize',10);
    hold off;
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ylabel('PSD exponent', 'FontSize', 24, 'FontWeight','normal');
    xlabel('{\itR} (z-transformed)', 'FontSize', 24, 'FontWeight','normal');
    title([task_name{task}, ' ', electrode,' - young'], 'FontSize', 24, 'FontWeight','normal');
    axis([1 5 0 2])
    [r, p] = corrcoef(R2_transf(abs(zscore(R2_transf)) <= 2.5),exponent_young(task+1, abs(zscore(R2_transf)) <= 2.5));
    xt = 3.5; yt = .5;
    str = {['r = ' num2str(r(1, 2))], ['p = ' num2str(p(1, 2))]};
    text(xt,yt,str,'FontSize',14)
    
    figure;
    plot(R2_transf(abs(zscore(R2_transf)) <= 2.5),offset_young(task+1, abs(zscore(R2_transf)) <= 2.5), 'o', 'color', 'k',  'MarkerSize',10); lsline
    hold on; plot(R2_transf(outliers),offset_young(task+1, outliers), 'o', 'color', 'k', 'MarkerFaceColor', 'r',  'MarkerSize',10);
    hold off;
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ylabel('PSD offset', 'FontSize', 24, 'FontWeight','normal');
    xlabel('{\itR} (z-transformed)', 'FontSize', 24, 'FontWeight','normal');
    title([task_name{task}, ' ', electrode,'  - young'], 'FontSize', 24, 'FontWeight','normal');
    axis([1 5 -.5 2]); yt = 0;
    [r, p] = corrcoef(R2_transf(abs(zscore(R2_transf)) <= 2.5),offset_young(task+1, abs(zscore(R2_transf)) <= 2.5));
    str = {['r = ' num2str(r(1, 2))], ['p = ' num2str(p(1, 2))]};
    text(xt,yt,str,'FontSize',16)
end

%% older
for task=1:2
    R2_transf = squeeze(fisherz(sqrt(r_squared_error_older(task+1, :, 1))));
    outliers = find(abs(zscore(R2_transf)) > 2.5);
    figure;
    plot(R2_transf(abs(zscore(R2_transf)) <= 2.5),exponent_older(task+1, abs(zscore(R2_transf)) <= 2.5), 'o', 'color', 'k',  'MarkerSize',10); lsline
    hold on; plot(R2_transf(outliers ),exponent_older(task+1, outliers ), 'o', 'color', 'k', 'MarkerFaceColor', 'r',  'MarkerSize',10);
    hold off;
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ylabel('Exponent', 'FontSize', 24, 'FontWeight','normal');
    xlabel('Fisher transformed R2', 'FontSize', 24, 'FontWeight','normal');
    title([task_name{task}, ' ', electrode,'  - older'], 'FontSize', 24, 'FontWeight','normal');
    axis([1 5 0 2])
    [r, p] = corrcoef(R2_transf(abs(zscore(R2_transf)) <= 2.5),exponent_older(task+1, abs(zscore(R2_transf)) <= 2.5));
    xt = 3.5; yt = .5;
    str = {['r = ' num2str(r(1, 2))], ['p = ' num2str(p(1, 2))]};
    text(xt,yt,str,'FontSize',16)
    
    figure;
    plot(R2_transf(abs(zscore(R2_transf)) <= 2.5),offset_older(task+1, abs(zscore(R2_transf)) <= 2.5), 'o', 'color', 'k',  'MarkerSize',10); lsline
    hold on; plot(R2_transf(outliers ),offset_older(task+1, outliers ), 'o', 'color', 'k', 'MarkerFaceColor', 'r',  'MarkerSize',10);
    hold off;
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ylabel('Offset', 'FontSize', 24, 'FontWeight','normal');
    xlabel('Fisher transformed R2', 'FontSize', 24, 'FontWeight','normal');
    title([task_name{task}, ' ', electrode,' - older'], 'FontSize', 24, 'FontWeight','normal');
    axis([1 5 -.5 2]); yt = 0;
    [r, p] = corrcoef(R2_transf(abs(zscore(R2_transf)) <= 2.5),offset_older(task+1, abs(zscore(R2_transf)) <= 2.5));
    str = {['r = ' num2str(r(1, 2))], ['p = ' num2str(p(1, 2))]};
    text(xt,yt,str,'FontSize',16)
end

%% both group together - scatter plots - R2 (transformed) vs exponent and offset
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
load r_squared_error_older % (task, p, chan, r2 and error)
load r_squared_error_young
load exponent_young; load exponent_older; load offset_young; load offset_older; 
chan = 16; % FCz = 16
task_name = {'Simple RT', 'Go/no-go'}; electrode = 'FCz';
for task=2 %1:2

    R2_transf_young = squeeze(fisherz(sqrt(r_squared_error_young(task, :, chan, 1)))); 
    R2_transf_older = squeeze(fisherz(sqrt(r_squared_error_older(task, :, chan, 1))));
    
%     % exclude outliers
%     incl_young = find(abs(zscore(R2_transf_young )) <= 2.5);
%     incl_older = find(abs(zscore(R2_transf_older )) <= 2.5);
    figure; box off
    plot([R2_transf_young; R2_transf_older],[squeeze(exponent_young(task, :, chan))'; squeeze(exponent_older(task, :, chan))'],...
        'o', 'color', 'w',  'MarkerSize',.1); h = lsline; hold on
    set(h(1),'color','k')
    plot(R2_transf_older,exponent_older(task, :, chan), ...
            'o', 'color', [1 .5 .5], ...
            'MarkerFaceColor',[1 .5 .5], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5); hold on    

    plot(R2_transf_young,exponent_young(task, :, chan), ...
        'o', 'color', [.8 .8 .8],'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5);
    hold off; box off
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ax.LineWidth = 2.5; 
    ylabel('PSD exponent', 'FontSize', 32, 'FontWeight','normal');
    xlabel('{\itR} (z-transformed)', 'FontSize', 28, 'FontWeight','normal');
%     title([task_name{task}, ' ', electrode], 'FontSize', 32, 'FontWeight','normal');
%     axis([1.5 5 0 2])
%     y = [exponent_young(task+1, abs(zscore(R2_transf(1:length(younger)))) <= 2.5)'; exponent_older(task+1, abs(zscore(R2_transf(length(younger)+1:end))) <= 2.5)'];
%     x = [ones(length(R2_transf(abs(zscore(R2_transf)) <= 2.5)), 1), R2_transf(abs(zscore(R2_transf)) <= 2.5)];
%     b = regress(y,x)  
%     
%     % for a difference in R2 between 0.99 and 0.95
%     x1 = 3.05; y1 = b(1) + b(2)*x1
%     x2 = 3.62; y2 = b(1) + b(2)*x2
%     y1-y2
%     x1_r2 = ifisherz(x1).^2
%     x2_r2 = ifisherz(x2).^2

    [r, p] = corrcoef([R2_transf_young; R2_transf_older],[squeeze(exponent_young(task, :, chan))'; squeeze(exponent_older(task, :, chan))']);
    xt = 3.5; yt = .5;
    str = {['\itr = ' num2str(r(1, 2),'%4.3f')], ['\itp = ' num2str(p(1, 2),'%4.3f')]};
    text(xt,yt,str,'FontSize',18)
    

% offset
    figure; box off
    plot([R2_transf_young; R2_transf_older],[squeeze(offset_young(task, :, chan))'; squeeze(offset_older(task, :, chan))'],...
        'o', 'color', 'w',  'MarkerSize',.1); h = lsline; hold on
    set(h(1),'color','k')
    plot(R2_transf_older,offset_older(task, :, chan), ...
        'o', 'color', [1 .5 .5],'MarkerFaceColor',[1 .5 .5], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5);  hold on    
    plot(R2_transf_young,offset_young(task, :, chan), ...
        'o', 'color', [.8 .8 .8],'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5);
    hold off; box off
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ax.LineWidth = 2.5; 
    ylabel('PSD offset', 'FontSize', 32, 'FontWeight','normal');
    xlabel('{\itR} (z-transformed)', 'FontSize', 28, 'FontWeight','normal');
%     title([task_name{task}, ' ', electrode], 'FontSize', 32, 'FontWeight','normal');
%     axis([1.5 5  0 2])
    
    [r, p] = corrcoef([R2_transf_young; R2_transf_older],[squeeze(offset_young(task, :, chan))'; squeeze(offset_older(task, :, chan))']);
    xt = 3.5; yt = .5;
    str = {['\itr = ' num2str(r(1, 2),'%4.3f')], ['\itp = ' num2str(p(1, 2),'%4.3f')]};
    text(xt,yt,str,'FontSize',18)
end


%% 
% find participants where peakw were not detected for each freq band
% for task = 1:3
%     for peak = 1:3 %theta alpha beta
%         no_peak_detected_young(task, peak) = length(find(oscillations_power_young(task, :, peak) == 0))/size(oscillations_power_young, 2);
%         no_peak_detected_older(task, peak) = length(find(oscillations_power_older(task, :, peak) == 0))/size(oscillations_power_older, 2);
%     end
% end


for task = 1:3
    for peak = 1:3 %theta alpha beta
        peak_detected_young(task, peak) = length(find(oscillations_power_young(task, :, peak) > 0))/size(oscillations_power_young, 2);
        peak_detected_older(task, peak) = length(find(oscillations_power_older(task, :, peak) > 0))/size(oscillations_power_older, 2);
    end
end

% how many participants had more than one peak in each freq band
% number_peaks_young = task x participant x freq band
for task = 1:3
    for peak = 1:3 %theta alpha beta
        peaks_more_than_two_young(task, peak) = length(find(number_peaks_young(task, :, peak) > 1))/size(oscillations_power_young, 2);
        peaks_more_than_two_older(task, peak) = length(find(number_peaks_older(task, :, peak) > 1))/size(oscillations_power_older, 2);
    end
end

%% how many participants had four peaks fitted
% number_peaks_young = task x participant x freq band
clear peaks_young
for task = 1:3
    for p = 1:size(number_peaks_young, 2)
        peaks_young(task, p) = sum(number_peaks_young(task, p, :), 3);
    end
end


% percentage of participants with four peaks per task
for task =1:3
    perc_peaks_young(task) = length(find( peaks_young(task, :) == 4))/size(number_peaks_young, 2);
end

%% older
clear peaks_older
for task = 1:3
    for p = 1:size(number_peaks_older, 2)
        peaks_older(task, p) = sum(number_peaks_older(task, p, :), 3);
    end
end

% percentage of participants with four peaks per task
for task =1:3
    perc_peaks_older(task) = length(find(peaks_older(task, :) == 4))/size(number_peaks_older, 2);
end


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
    figure; box off;
    % % plot a line at zero
    % plot([0 0],[0 11], '--', 'color', [0 0 0]);
    % hold on
    % detection data
%     loglog(x_axis, mean_data_young,'color', 'k', 'LineWidth',1.5)
    plot(x_axis, (mean_data_young),'color', 'k', 'LineWidth',1.5)
    hold on
    jbfill(x_axis',(mean_data_young+se_data_young), (mean_data_young-se_data_young), 'k','k', 0.1)
    % GNG data
    hold on
%     loglog(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
    plot(x_axis, (mean_data_older),'color', 'r', 'LineWidth',1.5)
    hold on
    jbfill(x_axis', (mean_data_older+se_data_older), (mean_data_older-se_data_older),'r','r', 0.1)
    hold off
    box off;
    ax = gca;
    ax.LineWidth = 2.5; 
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    axis([1 35 -inf inf])
    xlabel('Frequency (Hz)', 'FontSize', 36, 'FontWeight','normal')
    ylabel(y_axis_label, 'FontSize', 36, 'FontWeight','normal')
    title(graph_title, 'FontSize', 36, 'FontWeight','normal')

end


%% function to plot all data points 3 tasks 2 groups

function plot_all_data_points(data_grp1_task1, data_grp1_task2, data_grp1_task3, data_grp2_task1, data_grp2_task2, data_grp2_task3, y_label_text)

% young group
% define colour of data points
%     for i=1:length(data_grp1_task1)
%         cmap(i, :) = [floor(255/length(data_grp1_task1)*i)/255 floor(255/length(data_grp1_task1)*i)/255 floor(255/length(data_grp1_task1)*i)/255];
%     end

    % plot data for young group - simple RT and go/nogo task

        figure;
    for y=1:length(data_grp1_task1)
        plot([1 2 3]+rand*0.2-0.1, [data_grp1_task1(y) data_grp1_task2(y) data_grp1_task3(y)] ,'-o', 'color', [.8 .8 .8], ...
            'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end

    %plot the median line
    yMean=nanmean([data_grp1_task1' data_grp1_task2' data_grp1_task3'], 1);
    plot([1 2 3],[yMean(1) yMean(2) yMean(3)] ,'Color','k','LineWidth',1.5);
    plot([1-0.3 1+0.3],[yMean(1) yMean(1)] ,'Color','k','LineWidth',5);
    plot([2-0.3 2+0.3],[yMean(2) yMean(2)] ,'Color','k','LineWidth',5);
    plot([3-0.3 3+0.3],[yMean(3) yMean(3)] ,'Color','k','LineWidth',5);
    

    % older group
%     for i=1:length(data_grp2_task1)
%         cmap(i, :) = [floor(255/length(data_grp2_task1)*i)/255 0 0];
%     end
    % plot data for older group - simple RT and go/nogo task

    
%     for y=1:length(data_grp2_task1)
%         plot([4 5 6]+rand*0.2-0.1, [data_grp2_task1(y) data_grp2_task2(y) data_grp2_task3(y)] ,'-o', 'color', cmap(y,:), ...
%             'MarkerFaceColor',cmap(y,:), 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
%         hold on;  
%     end
    for y=1:length(data_grp2_task1)
        plot([4 5 6]+rand*0.2-0.1, [data_grp2_task1(y) data_grp2_task2(y) data_grp2_task3(y)] ,'-o', 'color', [1 .8 .8], ...
            'MarkerFaceColor',[1 .8 .8], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;  
    end

     %plot the median line
    yMean=nanmean([data_grp2_task1' data_grp2_task2' data_grp2_task3'], 1);
    plot([4 5 6],[yMean(1) yMean(2) yMean(3)] ,'Color','k','LineWidth',1.5);
    plot([4-0.3 4+0.3],[yMean(1) yMean(1)] ,'Color','k','LineWidth',5);
    plot([5-0.3 5+0.3],[yMean(2) yMean(2)] ,'Color','k','LineWidth',5);
    plot([6-0.3 6+0.3],[yMean(3) yMean(3)] ,'Color','k','LineWidth',5);
    

    % axes('XColor','none');
    hold off;
    axis([0 7 -inf inf]);
    ax = gca;
    c = ax.Color;
    ax.FontSize = 18;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel=[];
    ylabel(y_label_text, 'FontSize', 24, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end

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
    
    
    %% group 2
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
    ylabel(y_label_text, 'FontSize', 40, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end

%% function plot data 2 groups 1 task
function plot_all_data_2groups(data_grp1, data_grp2, y_label_text)
       figure;
    % plot data for group 1
        yMean=nanmean(data_grp1);
    y_se = nanstd(data_grp1)/sqrt(length(data_grp1));
    
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
     %plot the mean line    
    plot([1-0.3 1+0.3],[yMean yMean] ,'Color','k','LineWidth',5);
   

    %% group 2
    yMean=nanmean(data_grp2);
    y_se = nanstd(data_grp2)/sqrt(length(data_grp2));
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
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end


function plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)

    % plot data for young group - go/nogo task
    figure; box off; hold on
    
    % delete nan from data
    data_grp1_task1(isnan(data_grp1_task1)) = [];
    data_grp2_task1(isnan(data_grp2_task1)) = [];
    
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
    ax.LineWidth = 2.5; 
    ax.YAxis.FontSize = 24;
    ax.XAxis.FontSize = 32;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel= {'Young' 'Older'};
    xticks([1 2])
    if contains(y_label_text, 'z-transformed')
        ylabel(y_label_text, 'FontSize', 32, 'FontWeight','normal')
    else
        ylabel(y_label_text, 'FontSize', 40, 'FontWeight','normal')
    end
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end
