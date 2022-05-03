% regress out or adjust for spectral properties in passive pupil data out of
% variability in teh task-related pupil response - 26Nov2002 - Maria Ribeiro
clear; close all
% load pupil variability variables calculated in
% G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\pupil_std_percent_signal_change.m
load_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\';
load([load_dir 'std_pupil_avg_response']); % group x task
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
% plot_all_data_2tasks(std_pupil_avg_response{1,1}, std_pupil_avg_response{1,2},...
%     std_pupil_avg_response{2,1}, std_pupil_avg_response{2,2}, 'Pupil stdev')

% load PSD fooof variables for regression
% calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\foof_analysis_pupil_PSD_20s_long_epochs.m
filename = 'foof_results_pupil_20s_epochs.xlsx';
T = readtable([load_dir filename]);

% variable to save excel file with residuals
T_values = T.group;
variable_names = {'group', 'PDvar_simpleRT', 'PDvar_gng', 'PDvar_simpleRT_adjexp', 'PDvar_gng_adjexp', ...
    'PDvar_simpleRT_adjoffset', 'PDvar_gng_adjoffset', 'PDvar_simpleRT_adjslowfluct', 'PDvar_gng_adjslowfluct'};
T_values = [T_values, [std_pupil_avg_response{1, 1}; std_pupil_avg_response{2, 1}]]; % PD variability in simple RT task
T_values = [T_values, [std_pupil_avg_response{1, 2}; std_pupil_avg_response{2, 2}]]; % PD variability in gng task
task_name = {'Simple RT', 'Go/no-go'};
% Pupil variability adjusted for spectral exponent measured in passive task condition
Res = []; % std_pupil_avg_response = group x task
for task = 1:2
    Y = [std_pupil_avg_response{1, task}; std_pupil_avg_response{2, task}];
    X = [T.exponent_passive, ones(length(T.exponent_passive), 1)];
    [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
    % scatter plot and correlation analyses
    % scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt)
    scatter_and_corr(X(T.group==1), Y(T.group==1), 'PD SD', 'Exponent', ['Young - ', task_name{task}], 1, -inf)
    scatter_and_corr(X(T.group==2), Y(T.group==2), 'PD SD', 'Exponent', ['Older - ', task_name{task}], 1, -inf)
end


for task = 2 % go/no-go
    Y = [std_pupil_avg_response{1, task}; std_pupil_avg_response{2, task}];
    X = [T.exponent_passive, ones(length(T.exponent_passive), 1)];
    [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
    % scatter plot and correlation analyses
    % scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt)
    scatter_and_corr(X(T.group==1), Y(T.group==1), 'PD SD', 'PSD exponent', 'Young', 1, -inf)
    scatter_and_corr(X(T.group==2), Y(T.group==2), 'PD SD', 'PSD exponent', 'Older', 1, -inf)
end
%
%%plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(Res(T.group==1, 1), Res(T.group==1, 2),...
    Res(T.group==2, 1), Res(T.group==2, 2), 'PD variability residuals')

plot_all_data_onetask(Res(T.group==1, 2), Res(T.group==2, 2), 'PD SD residuals')
[H,P,CI,STATS] = ttest2(Res(T.group==1, 2), Res(T.group==2, 2))
T_values = [T_values, Res]; % exponent

%% test if exponent measured in simple RT or gng tasks has an effect
% Res = []; % std_pupil_avg_response = group x task
% for task = 1:2
%     Y = [std_pupil_avg_response{1, task}; std_pupil_avg_response{2, task}];
%     if task == 1
%         X = [T.exponent_simpleRT, ones(length(T.exponent_passive), 1)];
%     else
%         X = [T.exponent_gng, ones(length(T.exponent_passive), 1)];
%     end
%     [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
%     % scatter plot and correlation analyses
%     % scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt)
%     scatter_and_corr(X(T.group==1), Y(T.group==1), 'Pupil stdev', 'Exponent', ['Young - ', task_name{task}], 1, -inf)
%     scatter_and_corr(X(T.group==2), Y(T.group==2), 'Pupil stdev', 'Exponent', ['Older - ', task_name{task}], 1, -inf)
% end
% 
% % plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
% plot_all_data_2tasks(Res(T.group==1, 1), Res(T.group==1, 2),...
%     Res(T.group==2, 1), Res(T.group==2, 2), 'Pupil variability residuals')
%% Pupil variability adjusted for spectral offset
close all
Res = []; % std_pupil_avg_response = group x task
for task = 1:2
    Y = [std_pupil_avg_response{1, task}; std_pupil_avg_response{2, task}];
    X = [T.offset_passive, ones(length(T.offset_passive), 1)];
    [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
    % scatter plot and correlation analyses
    % scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt)
    scatter_and_corr(X(T.group==1), Y(T.group==1), 'PD SD', 'Offset', ['Young - ', task_name{task}], 1, -inf)
    scatter_and_corr(X(T.group==2), Y(T.group==2), 'PD SD', 'Offset', ['Older - ', task_name{task}], 1, -inf)
end

for task = 2 % go/no-go
    Y = [std_pupil_avg_response{1, task}; std_pupil_avg_response{2, task}];
    X = [T.offset_passive, ones(length(T.offset_passive), 1)];
    [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
    % scatter plot and correlation analyses
    % scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt)
    scatter_and_corr(X(T.group==1), Y(T.group==1), 'PD SD', 'PSD offset', 'Young', 1, -inf)
    scatter_and_corr(X(T.group==2), Y(T.group==2), 'PD SD', 'PSD offset', 'Older', 1, -inf)
end

% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(Res(T.group==1, 1), Res(T.group==1, 2),...
    Res(T.group==2, 1), Res(T.group==2, 2), 'PD variability residuals')

plot_all_data_onetask(Res(T.group==1, 2), Res(T.group==2, 2), 'PD SD residuals')
[H,P,CI,STATS] = ttest2(Res(T.group==1, 2), Res(T.group==2, 2))
T_values = [T_values, Res]; % offset

%%
% Res = []; % std_pupil_avg_response = group x task
% for task = 1:2
%     Y = [std_pupil_avg_response{1, task}; std_pupil_avg_response{2, task}];
%     if task == 1
%         X = [T.offset_simpleRT, ones(length(T.exponent_passive), 1)];
%     else
%         X = [T.offset_gng, ones(length(T.exponent_passive), 1)];
%     end
%     [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
%     % scatter plot and correlation analyses
%     % scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt)
%     scatter_and_corr(X(T.group==1), Y(T.group==1), 'Pupil stdev', 'Offset', ['Young - ', task_name{task}], 1, -inf)
%     scatter_and_corr(X(T.group==2), Y(T.group==2), 'Pupil stdev', 'Offset', ['Older - ', task_name{task}], 1, -inf)
% end
% 
% % plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
% plot_all_data_2tasks(Res(T.group==1, 1), Res(T.group==1, 2),...
%     Res(T.group==2, 1), Res(T.group==2, 2), 'Pupil variability residuals')
%% regressout power of slow fluctuations as measuerd in raw PSD
close all
% load data - data were calculated here:
% G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\pupil_spectrum_percent_signal_change.m
load([load_dir 'Pupil_PSD_20s_epochs_Young']); load([load_dir 'Pupil_PSD_20s_epochs_Older']); load([load_dir 'Frequencies_20s_epochs']);
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

% Res = [];
% for task = 1:2
%     Y = [std_pupil_avg_response{1, task}; std_pupil_avg_response{2, task}];
%     X = [[squeeze(10*log10(PSD_per_condition_young(task+1, :, 4)))'; squeeze(10*log10(PSD_per_condition_older(task+1, :, 4)))'],...
%         ones(size(PSD_per_condition_young, 2)+size(PSD_per_condition_older, 2), 1)];
%     [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
%     % scatter plot and correlation analyses
%     xlim_min = -inf;
%     scatter_and_corr(X(T.group==1), Y(T.group==1), 'Pupil stdev', '10xlog10 slow spectral power', ['Young - ', task_name{task}], 0, xlim_min)
%     scatter_and_corr(X(T.group==2), Y(T.group==2), 'Pupil stdev', '10xlog10 slow spectral power', ['Older - ', task_name{task}], 0, xlim_min)
% end
% 
% % plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
% plot_all_data_2tasks(Res(T.group==1, 1), Res(T.group==1, 2), Res(T.group==2, 1), Res(T.group==2, 2), 'Pupil variability residuals')
% 



% correlate with power of slow fluctuations in passive task
% frequency at average peak in corelation between raw passive PSD and pupil
% variability - avg peak = .4394 = freq to use [0.468750000000000]
% calculated in: G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability\corr_aperiodic_spect_Pupilstdev.m
Res = [];
for task = 2%1:2
    Y = [std_pupil_avg_response{1, task}; std_pupil_avg_response{2, task}];
    X = [[squeeze(log10(PSD_per_condition_young(1, :, 5)))'; squeeze(log10(PSD_per_condition_older(1, :, 5)))'],...
        ones(size(PSD_per_condition_young, 2)+size(PSD_per_condition_older, 2), 1)];
    [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
    % scatter plot and correlation analyses
    xlim_min = -inf;
    scatter_and_corr(X(T.group==1), Y(T.group==1), 'PD SD', 'log10(slow spectral power)', 'Young', 0, xlim_min)
    scatter_and_corr(X(T.group==2), Y(T.group==2), 'PD SD', 'log10(slow spectral power)', 'Older', 0, xlim_min)
end

% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(Res(T.group==1, 1), Res(T.group==1, 2), Res(T.group==2, 1), Res(T.group==2, 2), 'PD variability residuals')

plot_all_data_onetask(Res(T.group==1, 2), Res(T.group==2, 2), 'PD SD residuals')
[H,P,CI,STATS] = ttest2(Res(T.group==1, 2), Res(T.group==2, 2))
%%
T_values = [T_values, Res]; % slow PSD power
%
R = array2table(T_values, ...
    'VariableNames',variable_names);

filename = [load_dir filesep 'Pupilstdev_adjusted4exp_offset_slowpower.xlsx'];
writetable(R,filename,'Sheet',1,'Range','A1')

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
    
    % plot line at zero
    plot([0 6], [0 0], '--k');
    % axes('XColor','none');
    hold off;
    axis([0 6 -.035 .04]);
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

function scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt, xlim_min)
    [r, p] = corrcoef(X, Y);
    figure;
    plot(X, Y, 'o', 'color', 'k',  'MarkerSize',12, 'LineWidth', 1.5); 
    if p(1, 2)<.05
        lsline;
    end
    box off
    ax = gca; ax.FontSize = 24; ax.FontName = 'Arial';ax.Color = 'none';
    ax.LineWidth = 2.5; 
    ylabel(ylabel_txt, 'FontSize', 32, 'FontWeight','normal');
    xlabel(xlabel_txt, 'FontSize', 32, 'FontWeight','normal');
    title(title_txt, 'FontSize', 32, 'FontWeight','normal');
    axis([xlim_min inf -inf inf])
    xl = xlim; yl = ylim;
    xt = .1*(xl(2)-xl(1))/4+xl(1); yt = 3.5*(yl(2)-yl(1))/4+yl(1);
    str = {['r = ' num2str(r(1, 2))], ['p = ' num2str(p(1, 2))]}
    if incl_txt == 1
        text(xt,yt,str,'FontSize',14)
    end
end

function plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)

    % plot data for young group - go/nogo task
    figure; box off; hold on
    
    % plot line at zero
    plot([0 3], [0 0], '--k');
    
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
    axis([0 3 -.035 .04]);
    ax = gca;
    ax.LineWidth = 2.5; 
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 28;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel= {'Young' 'Older'};
    xticks([1 2])
    ylabel(y_label_text, 'FontSize', 28, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end
