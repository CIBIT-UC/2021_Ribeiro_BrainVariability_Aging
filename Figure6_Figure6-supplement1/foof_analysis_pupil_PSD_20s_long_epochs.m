% fooof analysis pupil
% onset
clear; close all
% load data - data were calculated here:
% G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\pupil_spectrum_long_epochs.m
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\');
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
settings.peak_width_limits = [-inf .25]; settings.peak_threshold = 3.5; settings.aperiodic_mode = 'fixed'; settings.min_peak_height = 0.2;
settings.max_n_peaks = 4; freq_range = [0.05 6];
% young
for task = 1:3
    for p = 1:size(PSD_per_condition_young, 2)
        % fooof_results = fooof(freqs, power_spectrum, f_range, settings,
        % return_model);[0.05 2] - frequency interval analysed
        fooof_results = fooof(Frequencies', squeeze(PSD_per_condition_young(task, p, :))', freq_range, settings, 1);
%         fooof_plot(fooof_results); title(['subj' num2str(p)]);
        exponent_young(task, p) = fooof_results.aperiodic_params(2);
        offset_young(task, p) = fooof_results.aperiodic_params(1);
        oscillations_young{task, p} = fooof_results.peak_params; % peaks detected: central frequency, power, band width
        r_squared_error_young(task, p, :) = [fooof_results.r_squared fooof_results.error];
        number_peaks_young(task, p) = size(fooof_results.peak_params, 1);
    end
end

% older
for task = 1:3
    for p = 1:size(PSD_per_condition_older, 2)
        count_theta = 0; count_alpha = 0; count_beta = 0;
        % fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
        fooof_results = fooof(Frequencies', squeeze(PSD_per_condition_older(task, p, :))', freq_range, settings, 1);
%         fooof_plot(fooof_results); title(['subj' num2str(p) ' - task ' num2str(task)]);
        exponent_older(task, p) = fooof_results.aperiodic_params(2);
        offset_older(task, p) = fooof_results.aperiodic_params(1);
        oscillations_older{task, p} = fooof_results.peak_params; % peaks detected

        r_squared_error_older(task, p, :) = [fooof_results.r_squared fooof_results.error];
        number_peaks_older(task, p) = size(fooof_results.peak_params, 1);
    end
end

%%stats

% [h,p,ci,stats] = ttest2(exponent_young(1, :), exponent_older(1, :))

%% peaks detected
% 
% % frequency bands defined based on grand averages
% % band1 - 0.05 - 0.25
% % band2 - 0.25 - 0.45
% % band3 - 0.45 - 0.75
% 
% 
% % create variable with frequency, power and bandwidth for each frequency
% % band for each participant, for each task
% 
% peaks_detected_young = zeros(size(oscillations_young, 2), 27);
% column_name = {'Band1_freq_pass' 'Band1_power_pass' 'Band1_width_pass'...
%     'Band2_freq_pass' 'Band2_power_pass' 'Band2_width_pass'...
%     'Band3_freq_pass' 'Band3_power_pass' 'Band3_width_pass'...
%     'Band1_freq_simpleRT' 'Band1_power_simpleRT' 'Band1_width_simpleRT'...
%     'Band2_freq_simpleRT' 'Band2_power_simpleRT' 'Band2_width_simpleRT'...
%     'Band3_freq_simpleRT' 'Band3_power_simpleRT' 'Band3_width_simpleRT'...
%         'Band1_freq_gng' 'Band1_power_gng' 'Band1_width_gng'...
%     'Band2_freq_gng' 'Band2_power_gng' 'Band2_width_gng'...
%     'Band3_freq_gng' 'Band3_power_gng' 'Band3_width_gng'};
% 
% % young group
% for p = 1:size(oscillations_young, 2)
%     for task = 1:3
%         if ~isempty(oscillations_young{task, p})
%             for pks = 1:size(oscillations_young{task, p}, 1)
%                 % use peak with higher power if two detected in same frequency band
%                 if oscillations_young{task, p}(pks, 1)<.25 && oscillations_young{task, p}(pks, 2) > peaks_detected_young(p, (task-1)*9+2) 
%                     peaks_detected_young(p, (task-1)*9+1:(task-1)*9+3) = oscillations_young{task, p}(pks, :);
%                 elseif oscillations_young{task, p}(pks, 1) >= .25 && oscillations_young{task, p}(pks, 1)< .45 && oscillations_young{task, p}(pks, 2) > peaks_detected_young(p, (task-1)*9+5) 
%                     peaks_detected_young(p, (task-1)*9+4:(task-1)*9+6) = oscillations_young{task, p}(pks, :);
%                 elseif oscillations_young{task, p}(pks, 1) >= .45 && oscillations_young{task, p}(pks, 2) > peaks_detected_young(p, (task-1)*9+8) 
%                     peaks_detected_young(p, (task-1)*9+7:(task-1)*9+9) = oscillations_young{task, p}(pks, :);
%                 end
%             end
%         end
%     end
% end
% 
% % older group
% peaks_detected_older = zeros(size(oscillations_older, 2), 27);
% for p = 1:size(oscillations_older, 2)
%     for task = 1:3
%         if ~isempty(oscillations_older{task, p})
%             for pks = 1:size(oscillations_older{task, p}, 1)
%                 % use peak with higher power if two detected in same frequency band
%                 if oscillations_older{task, p}(pks, 1)<.25 && oscillations_older{task, p}(pks, 2) > peaks_detected_older(p, (task-1)*9+2) 
%                     peaks_detected_older(p, (task-1)*9+1:(task-1)*9+3) = oscillations_older{task, p}(pks, :);
%                 elseif oscillations_older{task, p}(pks, 1) >= .25 && oscillations_older{task, p}(pks, 1)< .45 && oscillations_older{task, p}(pks, 2) > peaks_detected_older(p, (task-1)*9+5) 
%                     peaks_detected_older(p, (task-1)*9+4:(task-1)*9+6) = oscillations_older{task, p}(pks, :);
%                 elseif oscillations_older{task, p}(pks, 1) >= .45 && oscillations_older{task, p}(pks, 2) > peaks_detected_older(p, (task-1)*9+8) 
%                     peaks_detected_older(p, (task-1)*9+7:(task-1)*9+9) = oscillations_older{task, p}(pks, :);
%                 end
%             end
%         end
%     end
% end

%% create excel file - for SPSS analysis
T_values = [];
T_values =  [T_values, [ones(size(exponent_young, 2), 1); ones(size(exponent_older, 2), 1)*2]];
T_values =  [T_values, [exponent_young(1, :)'; exponent_older(1, :)']];
T_values =  [T_values, [exponent_young(2, :)'; exponent_older(2, :)']];
T_values =  [T_values, [exponent_young(3, :)'; exponent_older(3, :)']];

T_values =  [T_values, [offset_young(1, :)'; offset_older(1, :)']];
T_values =  [T_values, [offset_young(2, :)'; offset_older(2, :)']];
T_values =  [T_values, [offset_young(3, :)'; offset_older(3, :)']];

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

all_column_names = {'group','exponent_passive','exponent_simpleRT','exponent_gng', ...
    'offset_passive','offset_simpleRT','offset_gng', ...
    'R2_passive' 'R2_simpleRT', 'R2_gng',...
    'fisherR2Z_passive' 'fisherR2Z_simpleRT', 'fisherR2Z_gng',...
    'error_passive' 'error_simpleRT', 'error_gng'...
    'zscore_R2_passive' 'zscore_R2_simpleRT', 'zscore_R2_gng',...
    'zscore_fisherR2Z_passive' 'zscore_fisherR2Z_simpleRT', 'zscore_fisherR2Z_gng',...
    'zscore_error_passive' 'zscore_error_simpleRT', 'zscore_error_gng'};

T = array2table(T_values, ...
    'VariableNames',all_column_names);

filename = 'foof_results_pupil_20s_epochs.xlsx';
writetable(T,filename,'Sheet',1,'Range','A1')

% % find outliers
% outliers_young = find(abs(T.zscore_fisherR2Z_passive(T.group==1))>2.5);
% outliers_older = find(abs(T.zscore_fisherR2Z_passive(T.group==2))>2.5);
% 
% outliers_young = find(abs(T.zscore_fisherR2Z_simpleRT(T.group==1))>2.5);
% outliers_older = find(abs(T.zscore_fisherR2Z_simpleRT(T.group==2))>2.5);
% 
% outliers_young = find(abs(T.zscore_fisherR2Z_gng(T.group==1))>2.5);
% outliers_older = find(abs(T.zscore_fisherR2Z_gng(T.group==2))>2.5);
% % there is only one outlier in the gng in the older group!

%% plot_all_data_points(data_grp1_task1, data_grp1_task2, data_grp1_task3, data_grp2_task1, data_grp2_task2, data_grp2_task3, y_label_text)
load_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\';
% load PSD fooof variables for regression
% calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\foof_analysis_pupil_PSD_20s_long_epochs.m
filename = 'foof_results_pupil_20s_epochs.xlsx';
T = readtable([load_dir filename]);
exponent_young(1, :) = T.exponent_passive(T.group == 1);
exponent_older(1, :) = T.exponent_passive(T.group == 2);
offset_young(1, :) = T.offset_passive(T.group == 1);
offset_older(1, :) = T.offset_passive(T.group == 2);

r_squared_error_young(1, :, 1) = T.R2_passive(T.group == 1);
r_squared_error_older(1, :, 1) = T.R2_passive(T.group == 2);

r_squared_error_young(3, :, 1) = T.R2_gng(T.group == 1);
r_squared_error_older(3, :, 1) = T.R2_gng(T.group == 2);

plot_all_data_onetask(exponent_young(1, :), exponent_older(1, :), 'Pupil PSD exponent');
plot_all_data_onetask(offset_young(1, :), offset_older(1, :), 'Pupil PSD offset');

%% active tasks - simple RT vs gng
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(exponent_young(2, :), exponent_young(3, :), exponent_older(2, :), exponent_older(3, :), 'Exponent')

% plot_all_data_points(exponent_young(1, :), exponent_young(2, :), exponent_young(3, :),...
%     exponent_older(1, :), exponent_older(2, :), exponent_older(3, :), 'Exponent');

%% R2_square = (task, p, (R2, error))
plot_all_data_2groups(r_squared_error_young(1, :, 1), r_squared_error_older(1, :, 1),'R2');

plot_all_data_2tasks( r_squared_error_young(2, :, 1), r_squared_error_young(3, :, 1), r_squared_error_older(2, :, 1), r_squared_error_older(3, :, 1),  'R2')

plot_all_data_2tasks( fisherz(sqrt(r_squared_error_young(2, :, 1))), fisherz(sqrt(r_squared_error_young(3, :, 1))), fisherz(sqrt(r_squared_error_older(2, :, 1))), fisherz(sqrt(r_squared_error_older(3, :, 1))),  'FisherZ(R)')

%% passive task - used to calculate pupil ongoign fluctuations
plot_all_data_onetask(fisherz(sqrt(r_squared_error_young(1, :, 1))), fisherz(sqrt(r_squared_error_older(1, :, 1))), '{\itR} (z-transformed)');

[h,p,ci,stats] = ttest2(fisherz(sqrt(r_squared_error_young(1, :, 1))), fisherz(sqrt(r_squared_error_older(1, :, 1))))

%%
mean(r_squared_error_young(2, :, 1), 2)
mean(r_squared_error_young(3, :, 1), 2)
mean(r_squared_error_older(2, :, 1), 2)
mean(r_squared_error_older(3, :, 1), 2)


%% correlation between R2 and exponent and offset
task_name = {'Passive' 'Simple RT', 'Go/no-go'}; 
for task=1 

    R2_transf_young = squeeze(fisherz(sqrt(r_squared_error_young(task, :, 1)))); 
    R2_transf_older = squeeze(fisherz(sqrt(r_squared_error_older(task, :, 1))));
    
    % exclude outliers - there are no outliers!
    incl_young = find(abs(zscore(R2_transf_young )) <= 2.5);
    incl_older = find(abs(zscore(R2_transf_older )) <= 2.5);
    figure;
%     plot([R2_transf_young(incl_young); R2_transf_older(incl_older)],[squeeze(exponent_young(task, incl_young))'; squeeze(exponent_older(task, incl_older))'],...
%         'o', 'color', 'w',  'MarkerSize',.1); h = lsline; hold on
%     set(h(1),'color','k')
    plot(R2_transf_older(incl_older),exponent_older(task, incl_older), ...
        'o', 'color', [1 .5 .5],'MarkerFaceColor',[1 .5 .5], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5); hold on    
    plot(R2_transf_young(incl_young),exponent_young(task, incl_young), ...
        'o', 'color', [.8 .8 .8],'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5);
    hold off; box off
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ax.LineWidth = 2.5; 
    ylabel('PSD exponent', 'FontSize', 32, 'FontWeight','normal');
    xlabel('{\itR} (z-transformed)', 'FontSize', 28, 'FontWeight','normal');
%     title([task_name{task}], 'FontSize', 32, 'FontWeight','normal');
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
    
    [r, p] = corrcoef([R2_transf_young(incl_young); R2_transf_older(incl_older)],[squeeze(exponent_young(task, incl_young))'; squeeze(exponent_older(task, incl_older))']);
    xt = 3.5; yt = 1.95;
    str = {['\itr = ' num2str(r(1, 2),'%4.3f')], ['\itp = ' num2str(p(1, 2),'%4.3f')]};
    text(xt,yt,str,'FontSize',18)
    

% offset
    figure;
    plot([R2_transf_young(incl_young); R2_transf_older(incl_older)],[squeeze(offset_young(task, incl_young))'; squeeze(offset_older(task, incl_older))'],...
        'o', 'color', 'w',  'MarkerSize',.1); h = lsline; hold on
    set(h(1),'color','k')
    plot(R2_transf_older(incl_older),offset_older(task, incl_older), ...
         'o', 'color', [1 .5 .5],'MarkerFaceColor',[1 .5 .5], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5); hold on    
    plot(R2_transf_young(incl_young),offset_young(task, incl_young), ...
         'o', 'color', [.8 .8 .8],'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',12, 'LineWidth', 1.5);
    hold off; box off
    ax = gca; ax.FontSize = 18; ax.FontName = 'Arial';ax.Color = 'none';
    ax.LineWidth = 2.5; 
    ylabel('PSD offset', 'FontSize', 32, 'FontWeight','normal');
    xlabel('{\itR} (z-transformed)', 'FontSize', 28, 'FontWeight','normal');
%     title([task_name{task}], 'FontSize', 32, 'FontWeight','normal');
%     axis([1.5 5  0 2])
    
    [r, p] = corrcoef([R2_transf_young(incl_young); R2_transf_older(incl_older)],[squeeze(offset_young(task, incl_young))'; squeeze(offset_older(task, incl_older))']);
    xt = 3.3; yt = -4.2;
    str = {['\itr = ' num2str(r(1, 2),'%4.3f')], ['\itp = ' num2str(p(1, 2),'%4.3f')]};
    text(xt,yt,str,'FontSize',18)

end


%% error
% plot_all_data_points(r_squared_error_young(1, :, 2), r_squared_error_young(2, :, 2), r_squared_error_young(3, :, 2),...
%     r_squared_error_older(1, :, 2), r_squared_error_older(2, :, 2), r_squared_error_older(3, :, 2), 'Error');

plot_all_data_2groups(r_squared_error_young(1, :, 2), r_squared_error_older(1, :, 2),'Error');

plot_all_data_2tasks( r_squared_error_young(2, :, 2), r_squared_error_young(3, :, 2), r_squared_error_older(2, :, 2), r_squared_error_older(3, :, 2),  'Error')


%% 

% column_name = {'Band1_freq_pass' 'Band1_power_pass' 'Band1_width_pass'...
%     'Band2_freq_pass' 'Band2_power_pass' 'Band2_width_pass'...
%     'Band3_freq_pass' 'Band3_power_pass' 'Band3_width_pass'...
%     'Band1_freq_simpleRT' 'Band1_power_simpleRT' 'Band1_width_simpleRT'...
%     'Band2_freq_simpleRT' 'Band2_power_simpleRT' 'Band2_width_simpleRT'...
%     'Band3_freq_simpleRT' 'Band3_power_simpleRT' 'Band3_width_simpleRT'...
%         'Band1_freq_gng' 'Band1_power_gng' 'Band1_width_gng'...
%     'Band2_freq_gng' 'Band2_power_gng' 'Band2_width_gng'...
%     'Band3_freq_gng' 'Band3_power_gng' 'Band3_width_gng'};
% 
% 
% pupil_foof = readtable('foof_results_pupil_20s_epochs.xlsx');
% % for graphs change zeros to NaN
% 
% pupil_foof.Band1_power_pass(pupil_foof.Band1_power_pass == 0) = NaN;
% pupil_foof.Band2_power_pass(pupil_foof.Band2_power_pass == 0) = NaN;
% pupil_foof.Band3_power_pass(pupil_foof.Band3_power_pass == 0) = NaN;
% 
% pupil_foof.Band1_power_simpleRT(pupil_foof.Band1_power_simpleRT == 0) = NaN;
% pupil_foof.Band2_power_simpleRT(pupil_foof.Band2_power_simpleRT == 0) = NaN;
% pupil_foof.Band3_power_simpleRT(pupil_foof.Band3_power_simpleRT == 0) = NaN;
% 
% pupil_foof.Band1_power_gng(pupil_foof.Band1_power_gng == 0) = NaN;
% pupil_foof.Band2_power_gng(pupil_foof.Band2_power_gng == 0) = NaN;
% pupil_foof.Band3_power_gng(pupil_foof.Band3_power_gng == 0) = NaN;
% 
% 
% %% passive task
% % plot_all_data_2groups(pupil_foof.Band1_power_pass(pupil_foof.group == 1), pupil_foof.Band1_power_pass(pupil_foof.group == 2), 'Peak1 power')
% 
% plot_all_data_2groups(pupil_foof.Band2_power_pass(pupil_foof.group == 1), pupil_foof.Band2_power_pass(pupil_foof.group == 2), 'Peak2 power')
% 
% plot_all_data_2groups(pupil_foof.Band3_power_pass(pupil_foof.group == 1), pupil_foof.Band3_power_pass(pupil_foof.group == 2), 'Peak3 power')


%% simple RT vs gng
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
% plot_all_data_2tasks(pupil_foof.Band1_power_simpleRT(pupil_foof.group == 1), pupil_foof.Band1_power_gng(pupil_foof.group == 1),... 
%     pupil_foof.Band1_power_simpleRT(pupil_foof.group == 2), pupil_foof.Band1_power_gng(pupil_foof.group == 2), 'Peak 1 power')
% 
% plot_all_data_2tasks(pupil_foof.Band1_power_simpleRT(pupil_foof.group == 1), pupil_foof.Band1_power_gng(pupil_foof.group == 1),... 
%     pupil_foof.Band1_power_simpleRT(pupil_foof.group == 2), pupil_foof.Band1_power_gng(pupil_foof.group == 2), 'Peak 1 power')
% 
% plot_all_data_2tasks(pupil_foof.Band1_power_simpleRT(pupil_foof.group == 1), pupil_foof.Band1_power_gng(pupil_foof.group == 1),... 
%     pupil_foof.Band1_power_simpleRT(pupil_foof.group == 2), pupil_foof.Band1_power_gng(pupil_foof.group == 2), 'Peak 1 power')
%% plot_spectrum_mean_std_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)

title_text = {'Passive', 'Simple RT', 'Go/no-go'};
for task = 1:3
    plot_spectrum_mean_std_young_older(log10(squeeze(PSD_per_condition_young(task, :, :))), log10(squeeze(PSD_per_condition_older(task, :, :))), Frequencies, 6, title_text{task}, 'log10(PSD)')
end

%% plot graph function - group comparison
function plot_spectrum_mean_std_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)

    mean_data_young=mean(data_young(:, frequencies <= max_freq2plot & frequencies > 0), 1);
    se_data_young=std(data_young(:, frequencies <= max_freq2plot & frequencies > 0), 0, 1)/sqrt(size(data_young, 1));

    mean_data_older=mean(data_older(:, frequencies <= max_freq2plot & frequencies > 0), 1);
    se_data_older=std(data_older(:, frequencies <= max_freq2plot & frequencies > 0), 0, 1)/sqrt(size(data_older, 1));


    % color gng young = [0 .5 0]
    % color simple RT young =   [1 .5 0]
    % color gng older = [0 0 .5]
    % color simple RT older =  [.75 0.25 0]

    x_axis = frequencies(frequencies <= max_freq2plot & frequencies > 0); % frequencies
    figure;
    % % plot a line at zero
    % plot([0 0],[0 11], '--', 'color', [0 0 0]);
    % hold on
    % detection data
%     loglog(x_axis, mean_data_young,'color', 'k', 'LineWidth',1.5)
    plot(x_axis, (mean_data_young),'color', 'k', 'LineWidth',1.5)
    hold on
    jbfill(x_axis',(mean_data_young+se_data_young), (mean_data_young-se_data_young), 'k','k', 0.1)
    hold on
%     loglog(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
    plot(x_axis, (mean_data_older),'color', 'r', 'LineWidth',1.5)
    hold on
    jbfill(x_axis', (mean_data_older+se_data_older), (mean_data_older-se_data_older),'r','r', 0.1)
    hold off
    box off
    ax = gca; ax.LineWidth = 2.5; 
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
%     axis([1 35 -inf inf])
    xlabel('Frequency (Hz)', 'FontSize', 32, 'FontWeight','normal')
    ylabel(y_axis_label, 'FontSize', 32, 'FontWeight','normal')
%     title(graph_title, 'FontSize', 36, 'FontWeight','normal')

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
    ylabel(y_label_text, 'FontSize', 32, 'FontWeight','normal')
    
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
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel= [];
%     xticks([1 2])
    ylabel(y_label_text, 'FontSize', 32, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end


function plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)

    % plot data for young group - go/nogo task
    figure; box off; hold on
    
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
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.YAxis.FontSize = 24;
    ax.XTickLabel= {'Young' 'Older'};
    xticks([1 2])
    if contains(y_label_text, 'z-transformed')
        ylabel(y_label_text, 'FontSize', 32, 'FontWeight','normal')
        ax.XAxis.FontSize = 32;
    else
        ylabel(y_label_text, 'FontSize', 28, 'FontWeight','normal')
        ax.YAxis.FontSize = 18;
        ax.XAxis.FontSize = 28;
    end
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end
