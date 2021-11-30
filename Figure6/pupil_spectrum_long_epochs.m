clear; close all;
younger=[4	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46 50 51 53 54 56	59	62 66	68	72	74	76  78  80 81   82  84  85];
older=[7  8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60  61	63	64  65	67	69	70	71	73	75  77 79   83  86]; % AB60 was excluded due to noisy eeg - no alpha
% task={'W1', 'D1', 'D2', 'G1', 'G2'};
% centroparietal_cluster = [34 42 43 44 51]; %POz, P1, Pz, P2, CPz - for alpha calculation
group = {younger older};
std_pupil = cell(2, 5); std_pupil_avg_response = cell(2, 5);
% pupil_percent_signal_change = cell(2, 5); 
% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for grp = 1:2
    count = 0;
    for p = 8%group{grp}
        count = count + 1;
        pupil_directory=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\');
%         load([pupil_directory, 'RejectEpochsCorrect']); % data from G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\pupil_analysis_3b_RejectedEpochs_SaveRunData.m
        
        if intersect(p, [8, 11, 14, 15]) % in these participants use second repetition of recordign due to technical problems
            task = {'W2', 'D1', 'D2', 'G1', 'G2'};
        else
            task={'W1', 'D1', 'D2', 'G1', 'G2'};
        end
        
        clear pupil_data_percent_signal_change
        for t = 1:5 % task runs
        
            % calculate pupil go trials data
            filename=strcat(pupil_directory, 'AB', num2str(p), '_', task(t), '_PupilData.mat');
            load(filename{1})
            % use smoothed pupil diameter data with blink artefacts interpolated
            % measured in horizontal direction (x)
            
            % delete data at the beginning and end, before and after task finished 
            tmp = find(PupilData(:, 7) > 0); % trigger channel
            pupil_data_percent_signal_change = PupilData(tmp(1):tmp(end), 3)/mean(PupilData(tmp(1):tmp(end), 3)) - 1;
%            figure; plot(pupil_data_percent_signal_change); title(['AB' num2str(p)]); ylabel('Pupil %');
            % divide data into 20 sec epochs and calculate PSD for each
            % epoch - 20 sec = 20*240 time point
            for epoch = 1:floor(length(pupil_data_percent_signal_change)/(20*240))
                [Pxx(epoch, :), F] = pwelch(pupil_data_percent_signal_change((20*240)*(epoch-1)+1:(20*240)*epoch), [], [], [], 240);
                figure; plot(F(1:10), Pxx(epoch, 1:10)); title(['AB' num2str(p) 'epoch' num2str(epoch)]); ylabel('PSD');
%            figure; loglog(F(1:75), Pxx(1:75)); title(['AB' num2str(p)]); ylabel('PSD');
            end
            
            
            figure;
            plot(F(2:10), mean(Pxx(:, 2:10), 1)); title(['AB' num2str(p)]); ylabel('PSD');
            
            if grp == 1
                Pupil_PSD_20s_epochs_Young(t, count, :) =  mean(Pxx(:, 1:300), 1);
            else
                Pupil_PSD_20s_epochs_Older(t, count, :) = mean(Pxx(:, 1:300), 1);
            end
          
        end
    end
end

% save variables
Frequencies_20s_epochs = F(1:300);
save Frequencies_20s_epochs Frequencies_20s_epochs
save Pupil_PSD_20s_epochs_Young Pupil_PSD_20s_epochs_Young
save Pupil_PSD_20s_epochs_Older Pupil_PSD_20s_epochs_Older

%%
load Pupil_PSD_20s_epochs_Young
load Pupil_PSD_20s_epochs_Older
load Frequencies_20s_epochs

% plot data
run_name = {'Passive' 'Simple RT - 1' 'Simple RT - 2' 'Go/no-go - 1' 'Go/no-go - 2'};
for t = 1:5
    data_young = squeeze(Pupil_PSD_20s_epochs_Young(t, :, 2:end));
    data_older = squeeze(Pupil_PSD_20s_epochs_Older(t, :,  2:end));
    frequencies = Frequencies_20s_epochs(2:end); max_freq2plot = 2; graph_title = run_name(t); y_axis_label = 'Log10 Power (dB)';
    plot_spectrum_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)
    plot_loglog_spectrum_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)
end
%%
% task_name = {'Passive' 'Simple RT' 'Go/no-go'};
% 
% data_young_simple_RT = squeeze(mean(Pupil_PSD_Young(2:3, :, :)));
% data_young_gng = squeeze(mean(Pupil_PSD_Young(4:5, :, :)));
% data_older_simple_RT = squeeze(mean(Pupil_PSD_Older(2:3, :, :)));
% data_older_gng = squeeze(mean(Pupil_PSD_Older(4:5, :, :)));
% 
% frequencies = Frequencies; max_freq2plot = .75; graph_title = task_name(2); y_axis_label = 'Log10 Power (dB)';
% plot_spectrum_young_older(data_young_simple_RT, data_older_simple_RT, frequencies, max_freq2plot, graph_title, y_axis_label)
%  



%% plot graph function - group comparison
function plot_spectrum_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)

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
%     axis([1 35 -inf inf])
    xlabel('Frequency (Hz)', 'FontSize', 28, 'FontWeight','normal')
    ylabel(y_axis_label, 'FontSize', 28, 'FontWeight','normal')
    title(graph_title, 'FontSize', 28, 'FontWeight','normal')

end

function plot_loglog_spectrum_young_older(data_young, data_older, frequencies, max_freq2plot, graph_title, y_axis_label)

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
%     plot(x_axis, log10(mean_data_young),'color', 'k', 'LineWidth',1.5)
    hold on
%     jbfill(x_axis',log10(mean_data_young+se_data_young), log10(mean_data_young-se_data_young), 'k','k', 0.1)
    hold on
    loglog(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
%     plot(x_axis, log10(mean_data_older),'color', 'r', 'LineWidth',1.5)
    hold on
%     jbfill(x_axis', log10(mean_data_older+se_data_older), log10(mean_data_older-se_data_older),'r','r', 0.1)
    hold off
    ax = gca;
    % c = ax.Color;
    % legend('Detection', 'GNG')
    ax.FontSize = 20;
    ax.FontName = 'Arial';
    ax.Color = 'none';
%     axis([1 35 -inf inf])
    xlabel('Frequency (Hz)', 'FontSize', 28, 'FontWeight','normal')
    ylabel(y_axis_label, 'FontSize', 28, 'FontWeight','normal')
    title(graph_title, 'FontSize', 28, 'FontWeight','normal')

end
