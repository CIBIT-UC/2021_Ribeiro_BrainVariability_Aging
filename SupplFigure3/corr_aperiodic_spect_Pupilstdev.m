% correlations betwen amplitude of aperiod spectrum at each frequency and
% Pupil variability
clear; close all
% load pupil variability variables calculated in
% G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\pupil_std_percent_signal_change.m
load_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\';
load([load_dir 'std_pupil_avg_response']); % group x task

% load PSD fooof variables for regression 
% calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\foof_analysis_pupil_PSD_20s_long_epochs.m
filename = 'foof_results_pupil_20s_epochs.xlsx';
T = readtable([load_dir filename]);

% load data - data were calculated here:
% G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\pupil_spectrum_percent_signal_change.m
load Pupil_PSD_20s_epochs_Young; load Pupil_PSD_20s_epochs_Older; load Frequencies_20s_epochs;
Frequencies = Frequencies_20s_epochs;
% avg power spectrum across runs of same task condition
% group: 1 = young; 2 = older; task: 1 = passive; 2 = simple RT; 3 = gng
% Pupil_PSD_Young = task run x participant x frequencies
PSD_per_condition_young(1, :, :) = Pupil_PSD_20s_epochs_Young(1, :, :);
PSD_per_condition_young(2, :, :) = mean(Pupil_PSD_20s_epochs_Young(2:3, :, :), 1);
PSD_per_condition_young(3, :, :) = mean(Pupil_PSD_20s_epochs_Young(4:5, :, :), 1);

PSD_per_condition_older(1, :, :) = Pupil_PSD_20s_epochs_Older(1, :, :);
PSD_per_condition_older(2, :, :) = mean(Pupil_PSD_20s_epochs_Older(2:3, :, :), 1);
PSD_per_condition_older(3, :, :) = mean(Pupil_PSD_20s_epochs_Older(4:5, :, :), 1);

% frequencies fitted in fooof [0.05 2]
freqs = Frequencies_20s_epochs(Frequencies_20s_epochs>.05 & Frequencies_20s_epochs<2);
for grp = 1:2
    offset_passive = T.offset_passive(T.group == grp); exponent_passive = T.exponent_passive(T.group == grp);
    for p = 1:length(find(T.group == grp))
%         aperiodic_spectrum = fooof_results.aperiodic_params(1)-log10(k+fooof_results.freqs.^fooof_results.aperiodic_params(2));
        aperiodic_spectrum{grp}(p, :) = offset_passive(p)-log10(freqs.^exponent_passive(p));
    end
end

% [pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX,dataY,n_perm,tail,alpha_level,stat,reports,seed_state)
% by Groppe, D.M. - Each element of dataX should be paired with an element of dataY.
n_perm = 5000; tail = 0; alpha_level = .05; stat = 'linear'; reports = 0; 
group_name = {'Young', 'Older'}; task_name = {'Simple RT', 'Go/no-go'};
for grp = 1:2
    for task = 2%1:2
        % aperiodic fooof fitted spectrum
        dataX = aperiodic_spectrum{grp};
        dataY = repmat(std_pupil_avg_response{grp, task}, 1, size(dataX, 2));
        [pval(grp, task, :), corr_obs(grp, task, :), crit_corr(grp, task, :), ~, ~]=mult_comp_perm_corr(dataX,dataY);
        
        figure;ax = gca; ax.FontSize = 24; ax.FontName = 'Arial'; ax.Color = 'none';% c = ax.Color; 
        yyaxis left
        % plot grey background at freqs showing significant correlation
        sig_window = find(squeeze(corr_obs(grp, task, :)) > crit_corr(grp, task, 2))';
        for x = sig_window
            plot([freqs(x) freqs(x)],[0 .7], '-', 'color', [.8 .8 .8], 'LineWidth', 19 ); hold on
        end
        plot(freqs, squeeze(corr_obs(grp, task, :)), '-k', 'LineWidth', 2); hold on
%         % plot line at thrshold correlation r
%         plot(freqs, ones(length(freqs), 1)*crit_corr(grp, task, 2), '-k', 'LineWidth', .75); hold on
        
        set(gca,'ycolor','k') 
        xlabel('Frequency (Hz)', 'FontSize', 30)
        ylabel('Correlation \itr', 'FontSize', 30)
        axis([min(freqs), max(freqs), 0 .7])
        yyaxis right
        plot(freqs, squeeze(pval(grp, task, :)), '--', 'color', [0 158 115]/255, 'LineWidth', 2); hold on
        % plot line at p-value threshold
%         pvalue = squeeze(pval(grp, task, :));
%         plot([max(freqs(pvalue<=.05)) max(freqs(pvalue<.05))], [min(pvalue) max(pvalue)], '--r', 'LineWidth', .75);
        
        ylabel('P-value', 'FontSize', 30)
        set(gca,'ycolor', [0 158 115]/255) 
%         title([group_name{grp}, ' - ' ,task_name{task}], FontSize', 32, 'FontWeight','normal');
        title(group_name{grp}, 'FontSize', 32, 'FontWeight','normal');
        axis([min(freqs), max(freqs), -inf inf])
        
        % raw spectrum PSD in log10 scale
        if grp == 1
            dataX_PSD = squeeze(log10(PSD_per_condition_young(1, :, Frequencies_20s_epochs>0 & Frequencies_20s_epochs<2))); 
        else
            dataX_PSD = squeeze(log10(PSD_per_condition_older(1, :, Frequencies_20s_epochs>0 & Frequencies_20s_epochs<2))); 
        end
         dataY = repmat(std_pupil_avg_response{grp, task}, 1, size(dataX_PSD, 2));
        [pval_raw(grp, task, :), corr_obs_raw(grp, task, :), crit_corr_raw(grp, task, :), ~, ~]=mult_comp_perm_corr(dataX_PSD,dataY);
        
        x_freqs = Frequencies_20s_epochs(Frequencies_20s_epochs>0 & Frequencies_20s_epochs<2);
        figure;ax = gca; ax.FontSize = 24; ax.FontName = 'Arial'; c = ax.Color; ax.Color = 'none';
        yyaxis left
        % plot grey background at freqs showing significant correlation
        sig_window = find(squeeze(corr_obs_raw(grp, task, :)) > crit_corr_raw(grp, task, 2))';
        for x = sig_window
            plot([freqs(x) freqs(x)],[0 .8], '-', 'color', [.8 .8 .8], 'LineWidth', 19 ); hold on
        end
        plot(x_freqs, squeeze(corr_obs_raw(grp, task, :)), '-k', 'LineWidth', 2); hold on
%         % plot line at thrshold correlation r
%         plot(x_freqs, ones(length(x_freqs), 1)*crit_corr_raw(grp, task, 2), '-k', 'LineWidth', .75); hold on
        set(gca,'ycolor','k') 
        xlabel('Frequency (Hz)', 'FontSize', 30)
        ylabel('Correlation \itr', 'FontSize', 30)
        axis([min(freqs), max(freqs), 0 .7])
        yyaxis right
        plot(x_freqs, squeeze(pval_raw(grp, task, :)), '--', 'color', [0 158 115]/255, 'LineWidth', 2); hold on
        % plot line at p-value threshold
%         pvalue = squeeze(pval_raw(grp, task, :));
%         plot([max(freqs(pvalue<=.05)) max(freqs(pvalue<.05))], [min(pvalue) max(pvalue)], '--r', 'LineWidth', .75);
        
        ylabel('P-value', 'FontSize', 30)
        set(gca,'ycolor',[0 158 115]/255) 
%         title([group_name{grp}, ' - ' ,task_name{task}], 'FontSize', 32, 'FontWeight','normal');
        title(group_name{grp}, 'FontSize', 32, 'FontWeight','normal');
        axis([min(x_freqs), max(x_freqs), -inf inf])
    end
end

%%
for grp = 1:2
    for task = 1:2
        [M,I] = max(corr_obs(grp, task, :));
        freq_max(grp, task) = freqs(I)
        [M,I] = max(corr_obs_raw(grp, task, :));
        freq_max_raw(grp, task) = freqs(I)
    end
end
