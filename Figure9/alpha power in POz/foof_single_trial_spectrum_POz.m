
% fooof analysis electrode POz - pre-stimulus epoch 3.5 sec before cue
% onset
clear; close all
% participants id
younger=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];

save_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\single_trial_pre_stim_spectra_POz';

% load data - data were calculated here:
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\pre_stimulus_power_spectrum.m
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra';
load PowerSpectralDensity_Young; load PowerSpectralDensity_Older; load Frequencies;

% PowerSpectralDensity_Young = frequencies x trials
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs.mat');

% power spectrum at POz electrode
% find channel index
for c = 1:length(chanlocs)
    if strcmp(chanlocs(c).labels, 'POz')
        chan = c;
    end
end

%% avg power spectrum across trials
% % group: 1 = young; 2 = older; task: 1 = passive; 2 = simple RT; 3 = gng
% for p = 1:size(PowerSpectralDensity_Young, 1)
%     avg_power_spectrum_young(1, p, :) = mean(PowerSpectralDensity_Young{p, 1, chan}, 2);
%     avg_power_spectrum_young(2, p, :) = mean([PowerSpectralDensity_Young{p, 2, chan} PowerSpectralDensity_Young{p, 3, chan}], 2);
%     avg_power_spectrum_young(3, p, :) = mean([PowerSpectralDensity_Young{p, 4, chan} PowerSpectralDensity_Young{p, 5, chan}], 2);
% end
% 
% for p = 1:size(PowerSpectralDensity_Older, 1)
%     avg_power_spectrum_older(1, p, :) = mean(PowerSpectralDensity_Older{p, 1, chan}, 2);
%     avg_power_spectrum_older(2, p, :) =  mean([PowerSpectralDensity_Older{p, 2, chan} PowerSpectralDensity_Older{p, 3, chan}], 2);
%     avg_power_spectrum_older(3, p, :) = mean([PowerSpectralDensity_Older{p, 4, chan} PowerSpectralDensity_Older{p, 5, chan}], 2);
% end

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
settings.peak_width_limits = [1 inf]; settings.peak_threshold = 1; settings.aperiodic_mode = 'fixed'; settings.min_peak_height = 0.1;
settings.max_n_peaks = 2;
% young
for p = 1:size(PowerSpectralDensity_Young, 1)
    exponent = {}; oscillations_power = {}; oscillations_frequency = {}; oscillations_peak_width = {}; oscillations = {};
    r_squared_error = {}; number_peaks = {};
    for task = 1:5
        for trial = 1:size(PowerSpectralDensity_Young{p, task, chan}, 2)
            count_theta = 0; count_alpha = 0; count_beta = 0;
            % fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
            fooof_results = fooof(Frequencies', PowerSpectralDensity_Young{p, task, chan}(:, trial)', [1 35], settings, 1);
%             fooof_plot(fooof_results)
            exponent{task}(trial) = fooof_results.aperiodic_params(2);
            oscillations_power{task}(trial, :) = [0 0 0]; % theta alpha beta - if oscillations not detected power = 0;
            oscillations_frequency{task}(trial, :) = [NaN NaN NaN];
            oscillations_peak_width{task}(trial, :) = [NaN NaN NaN];
            for peak = 1:size(fooof_results.peak_params, 1)
                oscillations{p, task}(trial, peak, :) = fooof_results.peak_params(peak, :);
                if fooof_results.peak_params(peak, 1)>=4 && fooof_results.peak_params(peak, 1)< 8 % theta
                    count_theta = count_theta + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power{task}(trial, 1)
                        oscillations_power{task}(trial, 1) = fooof_results.peak_params(peak, 2); % power
                        oscillations_frequency{task}(trial, 1) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width{task}(trial, 1) = fooof_results.peak_params(peak, 3); % band width
                    end
                elseif fooof_results.peak_params(peak, 1)>=8 && fooof_results.peak_params(peak, 1)< 13 % alpha
                    count_alpha = count_alpha + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power{task}(trial, 2)
                        oscillations_power{task}(trial, 2) = fooof_results.peak_params(peak, 2);
                        oscillations_frequency{task}(trial, 2) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width{task}(trial, 2) = fooof_results.peak_params(peak, 3); % band width
                    end
                elseif fooof_results.peak_params(peak, 1)>=13 && fooof_results.peak_params(peak, 1)< 30 % beta
                    count_beta = count_beta + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power{task}(trial, 3)
                        oscillations_power{task}(trial, 3) = fooof_results.peak_params(peak, 2);
                        oscillations_frequency{task}(trial, 3) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width{task}(trial, 3) = fooof_results.peak_params(peak, 3); % band width
                    end
                end
            end
            r_squared_error{task}(trial, :) = [fooof_results.r_squared fooof_results.error];
            number_peaks{task}(trial, :) = [count_theta; count_alpha; count_beta];
        end
    end
    save([save_dir filesep 'young' filesep 'part_id_' num2str(younger(p))], 'exponent', 'oscillations_power', 'oscillations_frequency', 'oscillations_peak_width', 'oscillations',...
    'r_squared_error', 'number_peaks');
end

% older
for p = 1:size(PowerSpectralDensity_Older, 1)
    exponent = {}; oscillations_power = {}; oscillations_frequency = {}; oscillations_peak_width = {}; oscillations = {}
    r_squared_error = {}; number_peaks = {};
    for task = 1:5
        for trial = 1:size(PowerSpectralDensity_Older{p, task, chan}, 2)
            count_theta = 0; count_alpha = 0; count_beta = 0;
            % fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
            fooof_results = fooof(Frequencies', PowerSpectralDensity_Older{p, task, chan}(:, trial)', [1 35], settings, 1);
%             fooof_plot(fooof_results)
            exponent{task}(trial) = fooof_results.aperiodic_params(2);
            oscillations_power{task}(trial, :) = [0 0 0]; % theta alpha beta - if oscillations not detected power = 0;
            oscillations_frequency{task}(trial, :) = [NaN NaN NaN];
            oscillations_peak_width{task}(trial, :) = [NaN NaN NaN];
            for peak = 1:size(fooof_results.peak_params, 1)
                oscillations{p, task}(trial, peak, :) = fooof_results.peak_params(peak, :);
                if fooof_results.peak_params(peak, 1)>=4 && fooof_results.peak_params(peak, 1)< 8 % theta
                    count_theta = count_theta + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power{task}(trial, 1)
                        oscillations_power{task}(trial, 1) = fooof_results.peak_params(peak, 2); % power
                        oscillations_frequency{task}(trial, 1) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width{task}(trial, 1) = fooof_results.peak_params(peak, 3); % band width
                    end
                elseif fooof_results.peak_params(peak, 1)>=8 && fooof_results.peak_params(peak, 1)< 13 % alpha
                    count_alpha = count_alpha + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power{task}(trial, 2)
                        oscillations_power{task}(trial, 2) = fooof_results.peak_params(peak, 2);
                        oscillations_frequency{task}(trial, 2) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width{task}(trial, 2) = fooof_results.peak_params(peak, 3); % band width
                    end
                elseif fooof_results.peak_params(peak, 1)>=13 && fooof_results.peak_params(peak, 1)< 30 % beta
                    count_beta = count_beta + 1;
                    if fooof_results.peak_params(peak, 2) > oscillations_power{task}(trial, 3)
                        oscillations_power{task}(trial, 3) = fooof_results.peak_params(peak, 2);
                        oscillations_frequency{task}(trial, 3) = fooof_results.peak_params(peak, 1); % freq 
                        oscillations_peak_width{task}(trial, 3) = fooof_results.peak_params(peak, 3); % band width
                    end
                end
            end
            r_squared_error{task}(trial, :) = [fooof_results.r_squared fooof_results.error];
            number_peaks{task}(trial, :) = [count_theta; count_alpha; count_beta];
        end
    end
    save([save_dir filesep 'older' filesep 'part_id_' num2str(older(p))], 'exponent', 'oscillations_power', 'oscillations_frequency', 'oscillations_peak_width', 'oscillations',...
    'r_squared_error', 'number_peaks');
end

%% plot_all_data_points(data_grp1_task1, data_grp1_task2, data_grp1_task3, data_grp2_task1, data_grp2_task2, data_grp2_task3, y_label_text)
% exponent
% plot_all_data_2groups(data_grp1, data_grp2, y_label_text)
% passive task
plot_all_data_2groups(exponent{2}, exponent{3}, 'Exponent')

plot_all_data_2groups(oscillations_power{2}(:, 2), oscillations_power{3}(:, 2), 'Alpha power')



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
    plot(x_axis, log10(mean_data_young),'color', 'k', 'LineWidth',1.5)
    hold on
    jbfill(x_axis',log10(mean_data_young+se_data_young), log10(mean_data_young-se_data_young), 'k','k', 0.1)
    % GNG data
    hold on
%     loglog(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
    plot(x_axis, log10(mean_data_older),'color', 'r', 'LineWidth',1.5)
    hold on
    jbfill(x_axis', log10(mean_data_older+se_data_older), log10(mean_data_older-se_data_older),'r','r', 0.1)
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
%     plot([4 5],[yMean(1) yMean(2)] ,'Color','k','LineWidth',1.5);
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
