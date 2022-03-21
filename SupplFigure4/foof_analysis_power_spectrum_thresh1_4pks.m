% fooof analysis at chosen electrode - pre-stimulus epoch 3.5 sec before cue
% onset
clear; close all
electrode = input('Enter the name of electrode to analyse: ','s');

% load data - data were calculated here:
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\pre_stimulus_power_spectrum.m
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra';
load PowerSpectralDensity_Young; load PowerSpectralDensity_Older; load Frequencies;

load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs.mat');

% power spectrum at electrode
% find channel index
for c = 1:length(chanlocs)
    if strcmp(chanlocs(c).labels, electrode)
        chan = c;
    end
end

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
        exponent_young(task, p) = fooof_results.aperiodic_params(2);
        offset_young(task, p) = fooof_results.aperiodic_params(1);
        oscillations_power_young(task, p, :) = [0 0 0]; % theta alpha beta - if oscillations not detected power = 0;
        oscillations_frequency_young(task, p, :) = [NaN NaN NaN];
        oscillations_peak_width_young(task, p, :) = [NaN NaN NaN];
        for peak = 1:size(fooof_results.peak_params, 1)
            oscillations_young{p}(task, peak, :) = fooof_results.peak_params(peak, :);
            if fooof_results.peak_params(peak, 1)>=4 && fooof_results.peak_params(peak, 1)<=7  % theta - previous 8
                count_theta = count_theta + 1;
                % if there are two peaks within the same frequency band,use
                % peak with highest power
                if fooof_results.peak_params(peak, 2) > oscillations_power_young(task, p, 1)
                    oscillations_power_young(task, p, 1) = fooof_results.peak_params(peak, 2); % power
                    oscillations_frequency_young(task, p, 1) = fooof_results.peak_params(peak, 1); % freq 
                    oscillations_peak_width_young(task, p, 1) = fooof_results.peak_params(peak, 3); % band width
                end
            elseif fooof_results.peak_params(peak, 1)>7 && fooof_results.peak_params(peak, 1)< 14 % alpha previous 8 - 13
                count_alpha = count_alpha + 1;
                if fooof_results.peak_params(peak, 2) > oscillations_power_young(task, p, 2)
                    oscillations_power_young(task, p, 2) = fooof_results.peak_params(peak, 2);
                    oscillations_frequency_young(task, p, 2) = fooof_results.peak_params(peak, 1); % freq 
                    oscillations_peak_width_young(task, p, 2) = fooof_results.peak_params(peak, 3); % band width
                end
            elseif fooof_results.peak_params(peak, 1)>=14 && fooof_results.peak_params(peak, 1)< 30 % beta
                count_beta = count_beta + 1;
                if fooof_results.peak_params(peak, 2) > oscillations_power_young(task, p, 3)
                    oscillations_power_young(task, p, 3) = fooof_results.peak_params(peak, 2);
                    oscillations_frequency_young(task, p, 3) = fooof_results.peak_params(peak, 1); % freq 
                    oscillations_peak_width_young(task, p, 3) = fooof_results.peak_params(peak, 3); % band width
                end
            end
        end
        r_squared_error_young(task, p, :) = [fooof_results.r_squared fooof_results.error];
        number_peaks_young(task, p, :) = [count_theta; count_alpha; count_beta];
    end
end

% older
for task = 1:3
    for p = 1:size(avg_power_spectrum_older, 2)
        count_theta = 0; count_alpha = 0; count_beta = 0;
        % fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
        fooof_results = fooof(Frequencies', squeeze(avg_power_spectrum_older(task, p, :))', [1 35], settings, 1);
%         fooof_plot(fooof_results)
        exponent_older(task, p) = fooof_results.aperiodic_params(2);
        offset_older(task, p) = fooof_results.aperiodic_params(1);
        oscillations_power_older(task, p, :) = [0 0 0]; % theta alpha beta - if oscillations not detected power = 0;
        oscillations_frequency_older(task, p, :) = [NaN NaN NaN];
        oscillations_peak_width_older(task, p, :) = [NaN NaN NaN];
        for peak = 1:size(fooof_results.peak_params, 1)
            oscillations_older{p}(task, peak, :) = fooof_results.peak_params(peak, :);
            if fooof_results.peak_params(peak, 1)>=4 && fooof_results.peak_params(peak, 1)<= 7 % theta
                count_theta = count_theta + 1;
                if fooof_results.peak_params(peak, 2) > oscillations_power_older(task, p, 1)
                    oscillations_power_older(task, p, 1) = fooof_results.peak_params(peak, 2); % power
                    oscillations_frequency_older(task, p, 1) = fooof_results.peak_params(peak, 1); % freq 
                    oscillations_peak_width_older(task, p, 1) = fooof_results.peak_params(peak, 3); % band width
                end
            elseif fooof_results.peak_params(peak, 1)>7 && fooof_results.peak_params(peak, 1)< 14 % alpha
                count_alpha = count_alpha + 1;
                if fooof_results.peak_params(peak, 2) > oscillations_power_older(task, p, 2)
                    oscillations_power_older(task, p, 2) = fooof_results.peak_params(peak, 2);
                    oscillations_frequency_older(task, p, 2) = fooof_results.peak_params(peak, 1); % freq 
                    oscillations_peak_width_older(task, p, 2) = fooof_results.peak_params(peak, 3); % band width
                end
            elseif fooof_results.peak_params(peak, 1)>=14 && fooof_results.peak_params(peak, 1)< 30 % beta
                count_beta = count_beta + 1;
                if fooof_results.peak_params(peak, 2) > oscillations_power_older(task, p, 3)
                    oscillations_power_older(task, p, 3) = fooof_results.peak_params(peak, 2);
                    oscillations_frequency_older(task, p, 3) = fooof_results.peak_params(peak, 1); % freq 
                    oscillations_peak_width_older(task, p, 3) = fooof_results.peak_params(peak, 3); % band width
                end
            end
        end
        r_squared_error_older(task, p, :) = [fooof_results.r_squared fooof_results.error];
        number_peaks_older(task, p, :) = [count_theta; count_alpha; count_beta];
    end
end

%% adjusting exponent and offset for fitting quality - excluding outliers
clear exp_res offset_res incl_young incl_older
for task = 1:3
    R2_transf_young = squeeze(fisherz(sqrt(r_squared_error_young(task, :, 1))));
    incl_young{task} = find(abs(zscore(R2_transf_young)) <= 2.5); % without outliers
    R2_transf_older = squeeze(fisherz(sqrt(r_squared_error_older(task, :, 1))));
    incl_older{task} = find(abs(zscore(R2_transf_older)) <= 2.5);
    
    Y = [exponent_young(task, incl_young{task}), exponent_older(task, incl_older{task})]';
    X = [fisherz(sqrt([r_squared_error_young(task, incl_young{task}, 1), r_squared_error_older(task, incl_older{task}, 1)])), ones(length(Y), 1)];
    [B,BINT,exp_res{task},RINT,STATS] = regress(Y,X);
    
%     mdl = fitlm([exponent_young(task, incl_young{task}), exponent_older(task, incl_older{task})], fisherz(sqrt([r_squared_error_young(task, incl_young{task}, 1), r_squared_error_older(task, incl_older{task}, 1)])));
%     exp_res{task} = mdl.Residuals.Raw;

    Y = [offset_young(task, incl_young{task}), offset_older(task, incl_older{task})]';
    X = [fisherz(sqrt([r_squared_error_young(task, incl_young{task}, 1), r_squared_error_older(task, incl_older{task}, 1)])), ones(length(Y), 1)];
    [B,BINT,offset_res{task},RINT,STATS] = regress(Y,X);
    
%     mdl = fitlm([offset_young(task, incl_young{task}), offset_older(task, incl_older{task})], fisherz(sqrt([r_squared_error_young(task, incl_young{task}, 1), r_squared_error_older(task, incl_older{task}, 1)])));
%     offset_res{task} = mdl.Residuals.Raw;

end

% plot residuals
plot_all_data_2tasks(exp_res{2}(1:length(incl_young{2})), exp_res{3}(1:length(incl_young{3})), exp_res{2}(length(incl_young{2})+1:end), exp_res{3}(length(incl_young{3})+1:end), 'Exponent residuals')
plot_all_data_2tasks(offset_res{2}(1:length(incl_young{2})), offset_res{3}(1:length(incl_young{3})), offset_res{2}(length(incl_young{2})+1:end), offset_res{3}(length(incl_young{3})+1:end), 'Offset residuals')

%% add NAN in outliers in adjusted exponent and offset
exp_res_nan_young = nan(size(exponent_young, 2), 3);
exp_res_nan_older = nan(size(exponent_older, 2), 3);
offset_res_nan_young = nan(size(exponent_young, 2), 3);
offset_res_nan_older = nan(size(exponent_older, 2), 3);
for task = 1:3
     exp_res_nan_young(incl_young{task}, task) = exp_res{task}(1:length(incl_young{task}));
     exp_res_nan_older(incl_older{task}, task) = exp_res{task}(length(incl_young{task})+1:end);
     offset_res_nan_young(incl_young{task}, task) = offset_res{task}(1:length(incl_young{task}));
     offset_res_nan_older(incl_older{task}, task) = offset_res{task}(length(incl_young{task})+1:end);
end

%% create excel file - for SPSS analysis

T_values = [];
T_values(:, 1) = [ones(size(exponent_young, 2), 1); ones(size(exponent_older, 2), 1)*2];
T_values =  [T_values, [exponent_young(1, :)'; exponent_older(1, :)']];
T_values =  [T_values, [exponent_young(2, :)'; exponent_older(2, :)']];
T_values =  [T_values, [exponent_young(3, :)'; exponent_older(3, :)']];
T_values =  [T_values, [offset_young(1, :)'; offset_older(1, :)']];
T_values =  [T_values, [offset_young(2, :)'; offset_older(2, :)']];
T_values =  [T_values, [offset_young(3, :)'; offset_older(3, :)']];
% exponent and offset residuals after adjusting for model R2 and excluding
% outliers with bad fitting
T_values =  [T_values, [exp_res_nan_young(:, 1); exp_res_nan_older(:, 1)]];
T_values =  [T_values, [exp_res_nan_young(:, 2); exp_res_nan_older(:, 2)]];
T_values =  [T_values, [exp_res_nan_young(:, 3); exp_res_nan_older(:, 3)]];
T_values =  [T_values, [offset_res_nan_young(:, 1); offset_res_nan_older(:, 1)]];
T_values =  [T_values, [offset_res_nan_young(:, 2); offset_res_nan_older(:, 2)]];
T_values =  [T_values, [offset_res_nan_young(:, 3); offset_res_nan_older(:, 3)]];
% theta freq
T_values =  [T_values, [oscillations_frequency_young(1, :, 1)'; oscillations_frequency_older(1, :, 1)']];
T_values =  [T_values, [oscillations_frequency_young(2, :, 1)'; oscillations_frequency_older(2, :, 1)']];
T_values =  [T_values, [oscillations_frequency_young(3, :, 1)'; oscillations_frequency_older(3, :, 1)']];

%alpha_freq
T_values =  [T_values, [oscillations_frequency_young(1, :, 2)'; oscillations_frequency_older(1, :, 2)']];
T_values =  [T_values, [oscillations_frequency_young(2, :, 2)'; oscillations_frequency_older(2, :, 2)']];
T_values =  [T_values, [oscillations_frequency_young(3, :, 2)'; oscillations_frequency_older(3, :, 2)']];

%beta_freq
T_values =  [T_values, [oscillations_frequency_young(1, :, 3)'; oscillations_frequency_older(1, :, 3)']];
T_values =  [T_values, [oscillations_frequency_young(2, :, 3)'; oscillations_frequency_older(2, :, 3)']];
T_values =  [T_values, [oscillations_frequency_young(3, :, 3)'; oscillations_frequency_older(3, :, 3)']];

% theta power
T_values =  [T_values, [oscillations_power_young(1, :, 1)'; oscillations_power_older(1, :, 1)']];
T_values =  [T_values, [oscillations_power_young(2, :, 1)'; oscillations_power_older(2, :, 1)']];
T_values =  [T_values, [oscillations_power_young(3, :, 1)'; oscillations_power_older(3, :, 1)']];

%alpha_power
T_values =  [T_values, [oscillations_power_young(1, :, 2)'; oscillations_power_older(1, :, 2)']];
T_values =  [T_values, [oscillations_power_young(2, :, 2)'; oscillations_power_older(2, :, 2)']];
T_values =  [T_values, [oscillations_power_young(3, :, 2)'; oscillations_power_older(3, :, 2)']];

%beta_power
T_values =  [T_values, [oscillations_power_young(1, :, 3)'; oscillations_power_older(1, :, 3)']];
T_values =  [T_values, [oscillations_power_young(2, :, 3)'; oscillations_power_older(2, :, 3)']];
T_values =  [T_values, [oscillations_power_young(3, :, 3)'; oscillations_power_older(3, :, 3)']];

% theta band width
T_values =  [T_values, [oscillations_peak_width_young(1, :, 1)'; oscillations_peak_width_older(1, :, 1)']];
T_values =  [T_values, [oscillations_peak_width_young(2, :, 1)'; oscillations_peak_width_older(2, :, 1)']];
T_values =  [T_values, [oscillations_peak_width_young(3, :, 1)'; oscillations_peak_width_older(3, :, 1)']];

%alpha band width
T_values =  [T_values, [oscillations_peak_width_young(1, :, 2)'; oscillations_peak_width_older(1, :, 2)']];
T_values =  [T_values, [oscillations_peak_width_young(2, :, 2)'; oscillations_peak_width_older(2, :, 2)']];
T_values =  [T_values, [oscillations_peak_width_young(3, :, 2)'; oscillations_peak_width_older(3, :, 2)']];

%beta_peak_width
T_values =  [T_values, [oscillations_peak_width_young(1, :, 3)'; oscillations_peak_width_older(1, :, 3)']];
T_values =  [T_values, [oscillations_peak_width_young(2, :, 3)'; oscillations_peak_width_older(2, :, 3)']];
T_values =  [T_values, [oscillations_peak_width_young(3, :, 3)'; oscillations_peak_width_older(3, :, 3)']];

% number of peaks detected theta number_peaks_young(task, participant, peak)
T_values =  [T_values, [number_peaks_young(1, :, 1)'; number_peaks_older(1, :, 1)']];
T_values =  [T_values, [number_peaks_young(2, :, 1)'; number_peaks_older(2, :, 1)']];
T_values =  [T_values, [number_peaks_young(3, :, 1)'; number_peaks_older(3, :, 1)']];

% number of peaks detected alpha
T_values =  [T_values, [number_peaks_young(1, :, 2)'; number_peaks_older(1, :, 2)']];
T_values =  [T_values, [number_peaks_young(2, :, 2)'; number_peaks_older(2, :, 2)']];
T_values =  [T_values, [number_peaks_young(3, :, 2)'; number_peaks_older(3, :, 2)']];

% number of peaks detected beta
T_values =  [T_values, [number_peaks_young(1, :, 3)'; number_peaks_older(1, :, 3)']];
T_values =  [T_values, [number_peaks_young(2, :, 3)'; number_peaks_older(2, :, 3)']];
T_values =  [T_values, [number_peaks_young(3, :, 3)'; number_peaks_older(3, :, 3)']];

% model R2 
T_values =  [T_values, [r_squared_error_young(1, :, 1)'; r_squared_error_older(1, :, 1)']];
T_values =  [T_values, [r_squared_error_young(2, :, 1)'; r_squared_error_older(2, :, 1)']];
T_values =  [T_values, [r_squared_error_young(3, :, 1)'; r_squared_error_older(3, :, 1)']];

% model R2  - after fisher r to z
T_values =  [T_values, [fisherz(sqrt(r_squared_error_young(1, :, 1))); fisherz(sqrt(r_squared_error_older(1, :, 1)))]];
T_values =  [T_values, [fisherz(sqrt(r_squared_error_young(2, :, 1))); fisherz(sqrt(r_squared_error_older(2, :, 1)))]];
T_values =  [T_values, [fisherz(sqrt(r_squared_error_young(3, :, 1))); fisherz(sqrt(r_squared_error_older(3, :, 1)))]];

% model error
T_values =  [T_values, [r_squared_error_young(1, :, 2)'; r_squared_error_older(1, :, 2)']];
T_values =  [T_values, [r_squared_error_young(2, :, 2)'; r_squared_error_older(2, :, 2)']];
T_values =  [T_values, [r_squared_error_young(3, :, 2)'; r_squared_error_older(3, :, 2)']];

% model R2 - zscore
T_values =  [T_values, [zscore(r_squared_error_young(1, :, 1))'; zscore(r_squared_error_older(1, :, 1))']];
T_values =  [T_values, [zscore(r_squared_error_young(2, :, 1))'; zscore(r_squared_error_older(2, :, 1))']];
T_values =  [T_values, [zscore(r_squared_error_young(3, :, 1))'; zscore(r_squared_error_older(3, :, 1))']];

% model R2  - after fisher r to z - zscore
T_values =  [T_values, [zscore(fisherz(sqrt(r_squared_error_young(1, :, 1)))); zscore(fisherz(sqrt(r_squared_error_older(1, :, 1))))]];
T_values =  [T_values, [zscore(fisherz(sqrt(r_squared_error_young(2, :, 1)))); zscore(fisherz(sqrt(r_squared_error_older(2, :, 1))))]];
T_values =  [T_values, [zscore(fisherz(sqrt(r_squared_error_young(3, :, 1)))); zscore(fisherz(sqrt(r_squared_error_older(3, :, 1))))]];

% model error - zscore
T_values =  [T_values, [zscore(r_squared_error_young(1, :, 2))'; zscore(r_squared_error_older(1, :, 2))']];
T_values =  [T_values, [zscore(r_squared_error_young(2, :, 2))'; zscore(r_squared_error_older(2, :, 2))']];
T_values =  [T_values, [zscore(r_squared_error_young(3, :, 2))'; zscore(r_squared_error_older(3, :, 2))']];

% % variable with 1 if oscillations was detected 0 if it was not detected
% detected_peak_young  = oscillations_power_young;
% detected_peak_young(oscillations_power_young > 0) = 1;
% 
% detected_peak_older  = oscillations_power_older;
% detected_peak_older(oscillations_power_older > 0) = 1;
% 
% % theta detected
% T_values(:, 47) =  [detected_peak_young(1, :, 1)'; detected_peak_older(1, :, 1)'];
% T_values(:, 48) =  [detected_peak_young(2, :, 1)'; detected_peak_older(2, :, 1)'];
% T_values(:, 49) =  [detected_peak_young(3, :, 1)'; detected_peak_older(3, :, 1)'];
% 
% %alpha detected
% T_values(:, 50) =  [detected_peak_young(1, :, 2)'; detected_peak_older(1, :, 2)'];
% T_values(:, 51) =  [detected_peak_young(2, :, 2)'; detected_peak_older(2, :, 2)'];
% T_values(:, 52) =  [detected_peak_young(3, :, 2)'; detected_peak_older(3, :, 2)'];
% 
% %beta detected
% T_values(:, 53) =  [detected_peak_young(1, :, 3)'; detected_peak_older(1, :, 3)'];
% T_values(:, 54) =  [detected_peak_young(2, :, 3)'; detected_peak_older(2, :, 3)'];
% T_values(:, 55) =  [detected_peak_young(3, :, 3)'; detected_peak_older(3, :, 3)'];

T = array2table(T_values, ...
    'VariableNames',{'group','exponent_passive','exponent_simpleRT','exponent_gng',...
    'offset_passive','offset_simpleRT','offset_gng',...
    'exp_res_passive','exp_res_simpleRT','exp_res_gng',...
    'offset_res_passive','offset_res_simpleRT','offset_res_gng',...
    'theta_freq_passive' 'theta_freq_simpleRT', 'theta_freq_gng',...
    'alpha_freq_passive' 'alpha_freq_simpleRT', 'alpha_freq_gng',...
    'beta_freq_passive' 'beta_freq_simpleRT', 'beta_freq_gng',...
    'theta_power_passive' 'theta_power_simpleRT', 'theta_power_gng',...
    'alpha_power_passive' 'alpha_power_simpleRT', 'alpha_power_gng',...
    'beta_power_passive' 'beta_power_simpleRT', 'beta_power_gng',...
    'theta_peak_width_passive' 'theta_peak_width_simpleRT', 'theta_peak_width_gng',...
    'alpha_peak_width_passive' 'alpha_peak_width_simpleRT', 'alpha_peak_width_gng',...
    'beta_peak_width_passive' 'beta_peak_width_simpleRT', 'beta_peak_width_gng'...
    'theta_peaks_passive' 'theta_peaks_simpleRT', 'theta_peaks_gng',...
    'alpha_peaks_passive' 'alpha_peaks_simpleRT', 'alpha_peaks_gng',...
    'beta_peaks_passive' 'beta_peaks_simpleRT', 'beta_peaks_gng'...
    'R2_passive' 'R2_simpleRT', 'R2_gng',...
    'fisherR2Z_passive' 'fisherR2Z_simpleRT', 'fisherR2Z_gng',...
    'error_passive' 'error_simpleRT', 'error_gng'...
    'zscore_R2_passive' 'zscore_R2_simpleRT', 'zscore_R2_gng',...
    'zscore_fisherR2Z_passive' 'zscore_fisherR2Z_simpleRT', 'zscore_fisherR2Z_gng',...
    'zscore_error_passive' 'zscore_error_simpleRT', 'zscore_error_gng'...
    });


filename = ['foof_results_', electrode, '_thresh1_width8_4pks.xlsx'];
writetable(T,filename,'Sheet',1,'Range','A1')

%% read excel table to plot data
% T = readtable(filename);
% %%
% plot_all_data_2groups(T.R2_passive(T.group == 1 & T.R2_passive>.95), T.R2_passive(T.group == 2 & T.R2_passive>.95),'R2');
% %%
% plot_all_data_2tasks(T.R2_simpleRT(T.group == 1), T.R2_gng(T.group == 1),...
%     T.R2_simpleRT(T.group == 2), T.R2_gng(T.group == 2),  'R2')
% 
% plot_all_data_2groups(log(T.R2_passive(T.group == 1)), log(T.R2_passive(T.group == 2)),'R2');
% 
% plot_all_data_2tasks(log(T.R2_simpleRT(T.group == 1)),log(T.R2_gng(T.group == 1)),...
%     log(T.R2_simpleRT(T.group == 2)),log(T.R2_gng(T.group == 2)),  'R2')


%% plot_all_data_points(data_grp1_task1, data_grp1_task2, data_grp1_task3, data_grp2_task1, data_grp2_task2, data_grp2_task3, y_label_text)
% exponent
% plot_all_data_2groups(data_grp1, data_grp2, y_label_text)
% passive task
plot_all_data_2groups(exponent_young(1, :), exponent_older(1, :), 'Exponent ')
plot_all_data_2groups(offset_young(1, :), offset_older(1, :), 'Offset')
%%
plot_all_data_2groups(10.^exponent_young(1, :), 10.^exponent_older(1, :), 'Exponent')
plot_all_data_2groups(10.^offset_young(1, :), 10.^offset_older(1, :), 'Offset')

%% active tasks - simple RT vs gng
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(exponent_young(2, :), exponent_young(3, :), exponent_older(2, :), exponent_older(3, :), 'Exponent')
plot_all_data_2tasks(offset_young(2, :), offset_young(3, :), offset_older(2, :), offset_older(3, :), 'Offset')

% plot_all_data_points(exponent_young(1, :), exponent_young(2, :), exponent_young(3, :),...
%     exponent_older(1, :), exponent_older(2, :), exponent_older(3, :), 'Exponent');

%% R2_square
plot_all_data_2groups(r_squared_error_young(1, :, 1), r_squared_error_older(1, :, 1),'R2');

plot_all_data_2tasks( r_squared_error_young(2, :, 1), r_squared_error_young(3, :, 1), r_squared_error_older(2, :, 1), r_squared_error_older(3, :, 1),  'R2')

% plot_all_data_points(r_squared_error_young(1, :, 1), r_squared_error_young(2, :, 1), r_squared_error_young(3, :, 1),...
%     r_squared_error_older(1, :, 1), r_squared_error_older(2, :, 1), r_squared_error_older(3, :, 1), 'R2');

%% error
% plot_all_data_points(r_squared_error_young(1, :, 2), r_squared_error_young(2, :, 2), r_squared_error_young(3, :, 2),...
%     r_squared_error_older(1, :, 2), r_squared_error_older(2, :, 2), r_squared_error_older(3, :, 2), 'Error');

plot_all_data_2groups(r_squared_error_young(1, :, 2), r_squared_error_older(1, :, 2),'Error');

plot_all_data_2tasks( r_squared_error_young(2, :, 2), r_squared_error_young(3, :, 2), r_squared_error_older(2, :, 2), r_squared_error_older(3, :, 2),  'Error')

plot_all_data_2tasks( fisherz(sqrt(r_squared_error_young(2, :, 1))), fisherz(sqrt(r_squared_error_young(3, :, 1))), ...
    fisherz(sqrt(r_squared_error_older(2, :, 1))), fisherz(sqrt(r_squared_error_older(3, :, 1))),  'Fisher r to z')


%% scatter plots - R2 (transformed) vs exponent and offset
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
    ylabel('Exponent', 'FontSize', 24, 'FontWeight','normal');
    xlabel('Fisher transformed R2', 'FontSize', 24, 'FontWeight','normal');
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
    ylabel('Offset', 'FontSize', 24, 'FontWeight','normal');
    xlabel('Fisher transformed R2', 'FontSize', 24, 'FontWeight','normal');
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

%% both group together
task_name = {'Simple RT', 'Go/no-go'};
for task=1:2
    R2_transf = [squeeze(fisherz(sqrt(r_squared_error_young(task+1, :, 1)))); squeeze(fisherz(sqrt(r_squared_error_older(task+1, :, 1))))];
    outliers = find(abs(zscore(R2_transf)) > 2.5);
    figure;
    plot(R2_transf(abs(zscore(R2_transf)) <= 2.5),[exponent_young(task+1, abs(zscore(R2_transf(1:length(younger)))) <= 2.5)'; exponent_older(task+1, abs(zscore(R2_transf(length(younger)+1:end))) <= 2.5)'],...
        'o', 'color', 'k',  'MarkerSize',10); lsline
%     hold on; plot(R2_transf(outliers),exponent_young(task+1, outliers), 'o', 'color', 'k', 'MarkerFaceColor', 'r',  'MarkerSize',10);
    hold off;
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ylabel('Exponent', 'FontSize', 24, 'FontWeight','normal');
    xlabel('Fisher transformed R2', 'FontSize', 24, 'FontWeight','normal');
    title([task_name{task}, ' ', electrode,' - young'], 'FontSize', 24, 'FontWeight','normal');
    axis([1 5 0 2])
    y = [exponent_young(task+1, abs(zscore(R2_transf(1:length(younger)))) <= 2.5)'; exponent_older(task+1, abs(zscore(R2_transf(length(younger)+1:end))) <= 2.5)'];
    x = [ones(length(R2_transf(abs(zscore(R2_transf)) <= 2.5)), 1), R2_transf(abs(zscore(R2_transf)) <= 2.5)];
    b = regress(y,x)  
    
    % for a difference in R2 between 0.99 and 0.95
    x1 = 3.05; y1 = b(1) + b(2)*x1
    x2 = 3.62; y2 = b(1) + b(2)*x2
    y1-y2
    x1_r2 = ifisherz(x1).^2
    x2_r2 = ifisherz(x2).^2
    
    [r, p] = corrcoef(R2_transf(abs(zscore(R2_transf)) <= 2.5),exponent_young(task+1, abs(zscore(R2_transf)) <= 2.5));
    xt = 3.5; yt = .5;
    str = {['r = ' num2str(r(1, 2))], ['p = ' num2str(p(1, 2))]};
    text(xt,yt,str,'FontSize',14)
    
    figure;
    plot(R2_transf(abs(zscore(R2_transf)) <= 2.5),offset_young(task+1, abs(zscore(R2_transf)) <= 2.5), 'o', 'color', 'k',  'MarkerSize',10); lsline
    hold on; plot(R2_transf(outliers),offset_young(task+1, outliers), 'o', 'color', 'k', 'MarkerFaceColor', 'r',  'MarkerSize',10);
    hold off;
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ylabel('Offset', 'FontSize', 24, 'FontWeight','normal');
    xlabel('Fisher transformed R2', 'FontSize', 24, 'FontWeight','normal');
    title([task_name{task}, ' ', electrode,'  - young'], 'FontSize', 24, 'FontWeight','normal');
    axis([1 5 -.5 2]); yt = 0;
    [r, p] = corrcoef(R2_transf(abs(zscore(R2_transf)) <= 2.5),offset_young(task+1, abs(zscore(R2_transf)) <= 2.5));
    str = {['r = ' num2str(r(1, 2))], ['p = ' num2str(p(1, 2))]};
    text(xt,yt,str,'FontSize',16)
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


%%
% change zeros in power to NaN for graphs and SPSS stats
oscillations_power_young(oscillations_power_young == 0) = NaN;
oscillations_power_older(oscillations_power_older == 0) = NaN;
% alpha power
plot_all_data_points(oscillations_power_young(1, :, 2), oscillations_power_young(2, :, 2),oscillations_power_young(3, :, 2),...
   oscillations_power_older(1, :, 2), oscillations_power_older(2, :, 2), oscillations_power_older(3, :, 2), 'Alpha power');
%% plot power and frequency of alpha and beta oscillations
% plot_all_data_points(oscillations_power_young(1, :, 3), oscillations_power_young(2, :, 3),oscillations_power_young(3, :, 3),...
%    oscillations_power_older(1, :, 3), oscillations_power_older(2, :, 3), oscillations_power_older(3, :, 3), 'Beta power');

y_axis_text = {'Theta power' 'Alpha power', 'Beta power'};
y_axis_text_freq = {'Theta frequency (Hz)' 'Alpha frequency(Hz)', 'Beta frequency (Hz)'};

for peak = 2:3
    plot_all_data_2groups(squeeze(oscillations_power_young(1, :, peak)), squeeze(oscillations_power_older(1, :, peak)),y_axis_text{peak});
    plot_all_data_2tasks(squeeze(oscillations_power_young(2, :, peak )),squeeze(oscillations_power_young(3, :, peak )),...
        squeeze(oscillations_power_older(2, :, peak )), squeeze(oscillations_power_older(3, :, peak )),y_axis_text{peak})
    
    plot_all_data_2groups(squeeze(oscillations_frequency_young(1, :, peak)), squeeze(oscillations_frequency_older(1, :, peak)),y_axis_text_freq{peak});
    plot_all_data_2tasks(squeeze(oscillations_frequency_young(2, :, peak )),squeeze(oscillations_frequency_young(3, :, peak )),...
        squeeze(oscillations_frequency_older(2, :, peak )), squeeze(oscillations_frequency_older(3, :, peak )),y_axis_text_freq{peak})
end


%% plot_spectrum_mean_std_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)

title_text = {'Passive', 'Simple RT', 'Go/no-go'};
for task = 1:3
    plot_spectrum_mean_std_young_older(squeeze(avg_power_spectrum_young(task, :, :)), squeeze(avg_power_spectrum_older(task, :, :)), Frequencies, 35, title_text{task}, 'PSD (dB/Hz)')
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
    figure;
    % % plot a line at zero
    % plot([0 0],[0 11], '--', 'color', [0 0 0]);
    % hold on
    % detection data
%     loglog(x_axis, mean_data_young,'color', 'k', 'LineWidth',1.5)
    plot(x_axis, 10*log10(mean_data_young),'color', 'k', 'LineWidth',1.5)
    hold on
    jbfill(x_axis',10*log10(mean_data_young+se_data_young), 10*log10(mean_data_young-se_data_young), 'k','k', 0.1)
    % GNG data
    hold on
%     loglog(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
    plot(x_axis, 10*log10(mean_data_older),'color', 'r', 'LineWidth',1.5)
    hold on
    jbfill(x_axis', 10*log10(mean_data_older+se_data_older), 10*log10(mean_data_older-se_data_older),'r','r', 0.1)
    hold off
    ax = gca;
    % c = ax.Color;
    % legend('Detection', 'GNG')
    ax.FontSize = 20;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    axis([1 35 -inf inf])
    xlabel('Frequency (Hz)', 'FontSize', 28, 'FontWeight','normal')
    ylabel(y_axis_label, 'FontSize', 28, 'FontWeight','normal')
    title(graph_title, 'FontSize', 28, 'FontWeight','normal')

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
    ax.FontSize = 18;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel=[1 2 1 2];
    xticks([1 2 4 5])
    ylabel(y_label_text, 'FontSize', 24, 'FontWeight','normal')
    
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
