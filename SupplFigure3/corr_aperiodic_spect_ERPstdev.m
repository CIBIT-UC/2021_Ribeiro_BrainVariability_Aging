% correlations betwen amplitude of aperiod spectrum at each frequency and ERP
% variability
clear; close all
electrode = 'FCz'; %input('Enter the name of electrode to analyse: ','s');
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs.mat');
% power spectrum at electrode
% find channel index
for c = 1:length(chanlocs)
    if strcmp(chanlocs(c).labels, electrode)
        chan = c;
    end
end

% create variable with aperiodic spectrum from each participant
% load ERP variability variables calculated in
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\eeglab_analysis_5_channellevel_AllChannels_CorrectNoErr_StDev.m
load_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\';
load([load_dir 'ERP_avg_amp_stdev_young']);
load([load_dir 'ERP_avg_amp_stdev_older']);

% load PSD fooof variables for regression - FCz fitting with 4 peaks
% calculated in foof_analysis_power_spectrum_thresh1_4pks.m
table_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\';
load([table_dir, 'Frequencies.mat']); % load frequencies used in fooof analysis
load([table_dir, 'PowerSpectralDensity_Young']); load([table_dir, 'PowerSpectralDensity_Older']);
freqs = Frequencies(Frequencies>1 & Frequencies<35);
filename = ['foof_results_', electrode, '_thresh1_width8_4pks.xlsx'];
T = readtable([table_dir filename]);
% % exclude outliers with R2 z-score (calculated in the fisher r-to-z
% % transformed) >2.5
% outliers_young = find(abs(T.zscore_fisherR2Z_simpleRT(T.group==1))>2.5);
% outliers_older = find(abs(T.zscore_fisherR2Z_simpleRT(T.group==2))>2.5);
% T(abs(T.zscore_fisherR2Z_simpleRT)>2.5, :) = [];
% ERP_avg_amp_stdev_young(:, outliers_young, :) = [];
% ERP_avg_amp_stdev_older(:, outliers_older, :) = [];

% avg power spectrum across trials
% group: 1 = young; 2 = older; task: 1  = simple RT; 2 = gng
for p = 1:size(PowerSpectralDensity_Young, 1)
    avg_power_spectrum_young(1, p, :) = mean([PowerSpectralDensity_Young{p, 2, chan} PowerSpectralDensity_Young{p, 3, chan}], 2);
    avg_power_spectrum_young(2, p, :) = mean([PowerSpectralDensity_Young{p, 4, chan} PowerSpectralDensity_Young{p, 5, chan}], 2);
end

for p = 1:size(PowerSpectralDensity_Older, 1)
    avg_power_spectrum_older(1, p, :) =  mean([PowerSpectralDensity_Older{p, 2, chan} PowerSpectralDensity_Older{p, 3, chan}], 2);
    avg_power_spectrum_older(2, p, :) = mean([PowerSpectralDensity_Older{p, 4, chan} PowerSpectralDensity_Older{p, 5, chan}], 2);
end

% % remove outliers
% avg_power_spectrum_young(:, outliers_young, :) = [];
% avg_power_spectrum_older(:, outliers_older, :) = [];

for grp = 1:2
    offset_simpleRT = T.offset_simpleRT(T.group == grp); exponent_simpleRT = T.exponent_simpleRT(T.group == grp);
    offset_gng = T.offset_gng(T.group == grp); exponent_gng = T.exponent_gng(T.group == grp);
    for p = 1:length(find(T.group == grp))
%         aperiodic_spectrum = fooof_results.aperiodic_params(1)-log10(k+fooof_results.freqs.^fooof_results.aperiodic_params(2));
        aperiodic_spectrum{grp}(1, p, :) = offset_simpleRT(p)-log10(freqs.^exponent_simpleRT(p));
        aperiodic_spectrum{grp}(2, p, :) = offset_gng(p)-log10(freqs.^exponent_gng(p));
    end
end

%% [pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX,dataY,n_perm,tail,alpha_level,stat,reports,seed_state)
% by Groppe, D.M. - Each element of dataX should be paired with an element of dataY.
n_perm = 5000; tail = 0; alpha_level = .05; stat = 'linear'; reports = 0; 
group_name = {'Young', 'Older'}; task_name = {'Simple RT', 'Go/no-go'};
for grp = 1:2
    for task = 2%1:2
        % aperiodic fooof fitted spectrum
        dataX = squeeze(aperiodic_spectrum{grp}(task, :, :));
        if grp == 1
            dataY = repmat(squeeze(ERP_avg_amp_stdev_young(task, :, chan))', 1, size(dataX, 2));
        else
            dataY = repmat(squeeze(ERP_avg_amp_stdev_older(task, :, chan))', 1, size(dataX, 2));
        end
        [pval(grp, task, :), corr_obs(grp, task, :), crit_corr(grp, task, :), ~, ~]=mult_comp_perm_corr(dataX,dataY);
        
        figure;ax = gca; ax.FontSize = 24; ax.FontName = 'Arial'; c = ax.Color; ax.Color = 'none';
        yyaxis left
        % plot grey background at freqs showing significant correlation
        sig_window = find(squeeze(corr_obs(grp, task, :)) > crit_corr(grp, task, 2))';
        for x = sig_window
            plot([freqs(x) freqs(x)],[0 .95], '-', 'color', [.8 .8 .8], 'LineWidth', 10 ); hold on
        end
        plot(freqs, squeeze(corr_obs(grp, task, :)), '-k', 'LineWidth', 2); hold on
%         % plot line at thrshold correlation r
%         plot(freqs, ones(length(freqs), 1)*crit_corr(grp, task, 2), '-k', 'LineWidth', .75); hold on
        set(gca,'ycolor','k');
        xlabel('Frequency (Hz)', 'FontSize', 30)
        ylabel('Correlation \itr', 'FontSize', 30)
        axis([min(freqs), max(freqs), 0 .95])
        yyaxis right
        plot(freqs, squeeze(pval(grp, task, :)), '--', 'color',[0 158 115]/255, 'LineWidth', 2); hold on
        % plot line at p-value threshold
%         pvalue = squeeze(pval(grp, task, :));
%         plot([max(freqs(pvalue<=.05)) max(freqs(pvalue<.05))], [min(pvalue) max(pvalue)], '--r', 'LineWidth', .75);
        
        ylabel('P-value', 'FontSize', 30)
        set(gca,'ycolor',[0 158 115]/255) 
%         title([group_name{grp}, ' - ' ,task_name{task}], 'FontSize', 32, 'FontWeight','normal');
        title(group_name{grp}, 'FontSize', 32, 'FontWeight','normal');
        axis([min(freqs), max(freqs), -inf inf])
        
        % raw spectrum PSD in log10 scale
        if grp == 1
            dataX_PSD = squeeze(log10(avg_power_spectrum_young(task, :, Frequencies>1 & Frequencies<35))); 
        else
            dataX_PSD = squeeze(log10(avg_power_spectrum_older(task, :, Frequencies>1 & Frequencies<35))); 
        end
        [pval_raw(grp, task, :), corr_obs_raw(grp, task, :), crit_corr_raw(grp, task, :), ~, ~]=mult_comp_perm_corr(dataX_PSD,dataY);
        
        figure;ax = gca; ax.FontSize = 24; ax.FontName = 'Arial'; c = ax.Color; ax.Color = 'none';
        yyaxis left
        % plot grey background at freqs showing significant correlation
        sig_window = find(squeeze(corr_obs_raw(grp, task, :)) > crit_corr_raw(grp, task, 2))';
        for x = sig_window
            plot([freqs(x) freqs(x)],[0 .95], '-', 'color', [.8 .8 .8], 'LineWidth', 10 ); hold on
        end
        plot(freqs, squeeze(corr_obs_raw(grp, task, :)), '-k', 'LineWidth', 2); hold on
%         % plot line at threshold correlation r
%         plot(freqs, ones(length(freqs), 1)*crit_corr_raw(grp, task, 2), '-k', 'LineWidth', .75); hold on
        set(gca,'ycolor','k') 
        xlabel('Frequency (Hz)', 'FontSize', 30)
        ylabel('Correlation \itr', 'FontSize', 30)
        axis([min(freqs), max(freqs), 0 .95])
        yyaxis right
        plot(freqs, squeeze(pval_raw(grp, task, :)), '--', 'color', [0 158 115]/255, 'LineWidth', 2); hold on
        % plot line at p-value threshold
%         pvalue = squeeze(pval_raw(grp, task, :));
%         plot([max(freqs(pvalue<=.05)) max(freqs(pvalue<.05))], [min(pvalue) max(pvalue)], '--r', 'LineWidth', .75);
        
        ylabel('P-value', 'FontSize', 30)
        set(gca,'ycolor',[0 158 115]/255) 
%         title([group_name{grp}, ' - ' ,task_name{task}], 'FontSize', 32, 'FontWeight','normal');
        title(group_name{grp}, 'FontSize', 32, 'FontWeight','normal');
        axis([min(freqs), max(freqs), -inf inf])
    end
end

