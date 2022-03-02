% regress out or adjust for spectral properties in pre-stim EEG out of ERP
% variability - 26Nov2002 - Maria Ribeiro
clear; close all
% load ERP variability variables calculated in
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\eeglab_analysis_5_channellevel_AllChannels_CorrectNoErr_StDev.m

load_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\';

load([load_dir 'ERP_avg_amp_stdev_young']);
load([load_dir 'ERP_avg_amp_stdev_older']);


% load PSD fooof variables for regression - FCz fitting with 4 peaks
% calculated in foof_analysis_power_spectrum_FCz_thresh1_4pks
table_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\';
filename = 'foof_results_FCz_thresh1_width8_4pks.xlsx';
T = readtable([table_dir filename]);
% exclude outliers with R2 z-score (calculated in the fisher r-to-z
% transformed) >2.5
outliers_young = find(abs(T.zscore_fisherR2Z_simpleRT(T.group==1))>2.5);
outliers_older = find(abs(T.zscore_fisherR2Z_simpleRT(T.group==2))>2.5);
T(abs(T.zscore_fisherR2Z_simpleRT)>2.5, :) = [];
ERP_avg_amp_stdev_young(:, outliers_young, :) = [];
ERP_avg_amp_stdev_older(:, outliers_older, :) = [];

% variable to save excel file with residuals
T_values = T.group;
variable_names = {'group','exponent_simpleRT','exponent_gng',...
    'offset_simpleRT','offset_gng',...
    'slowpower_simpleRT','slowpower_gng'};
% ERP variability adjusted for spectral exponent - FCz - channel = 16
Res = [];
task = 1;
Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
X = [T.exponent_simpleRT, ones(length(T.exponent_simpleRT), 1)];
[B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
%% scatter plot and correlation analyses
% scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt)
scatter_and_corr(X(T.group==1, 1), Y(T.group==1), 'CNV SD', 'Exponent', 'Young - simple RT', 1, -inf)
scatter_and_corr(X(T.group==2, 1), Y(T.group==2), 'CNV SD', 'Exponent', 'Older - simple RT', 1, -inf)
%% go/no-go
task = 2;
Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
X = [T.exponent_gng, ones(length(T.exponent_gng), 1)];
[B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
% scatter plot and correlation analyses
scatter_and_corr(X(T.group==1), Y(T.group==1), 'CNV SD', 'PSD exponent', 'Young', 1, -inf)
scatter_and_corr(X(T.group==2), Y(T.group==2), 'CNV SD', 'PSD exponent', 'Older', 1, -inf)
%%
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(Res(1:size(ERP_avg_amp_stdev_young, 2), 1), Res(1:size(ERP_avg_amp_stdev_young, 2), 2), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 1), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'CNV variability residuals')
% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(Res(1:size(ERP_avg_amp_stdev_young, 2), 2),  Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'CNV SD residuals')

[H,P,CI,STATS] = ttest2(Res(1:size(ERP_avg_amp_stdev_young, 2), 2),  Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2))

T_values = [T_values, Res];

%% ERP variability adjusted for spectral offset - FCz - channel = 16
Res = [];
task = 1;
Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
X = [T.offset_simpleRT, ones(length(T.offset_simpleRT), 1)];
[B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
% scatter plot and correlation analyses
scatter_and_corr(X(T.group==1), Y(T.group==1), 'CNV SD', 'PSD offset', 'Young - simple RT', 0, -inf)
scatter_and_corr(X(T.group==2), Y(T.group==2), 'CNV SD', 'PSD offset', 'Older - simple RT', 0, -inf)
%% go/no-go
task = 2;
Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
X = [T.offset_gng, ones(length(T.offset_gng), 1)];
[B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
% scatter plot and correlation analyses
scatter_and_corr(X(T.group==1), Y(T.group==1), 'CNV SD', 'PSD offset', 'Young', 1, -inf)
scatter_and_corr(X(T.group==2), Y(T.group==2), 'CNV SD', 'PSD offset', 'Older', 1, -inf)
%%
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(Res(1:size(ERP_avg_amp_stdev_young, 2), 1), Res(1:size(ERP_avg_amp_stdev_young, 2), 2), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 1), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'CNV variability residuals')
% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(Res(1:size(ERP_avg_amp_stdev_young, 2), 2),  Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'CNV SD residuals')
[H,P,CI,STATS] = ttest2(Res(1:size(ERP_avg_amp_stdev_young, 2), 2),  Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2))
T_values = [T_values, Res];
%% regressout alpha power from ERP variability of both groups
Res = [];
task = 1;
Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
X = [T.alpha_power_simpleRT, ones(length(T.alpha_power_simpleRT), 1)];
[B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
% scatter plot and correlation analyses
xlim_min = -.2;
scatter_and_corr(X(T.group==1), Y(T.group==1), 'CNV SD', 'Alpha power', 'Young - simple RT', 0, xlim_min)
scatter_and_corr(X(T.group==2), Y(T.group==2), 'CNV SD', 'Alpha power', 'Older - simple RT', 0, xlim_min)
task = 2;
Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
X = [T.alpha_power_gng, ones(length(T.alpha_power_simpleRT), 1)];
[B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
% scatter plot and correlation analyses go/no-go
scatter_and_corr(X(T.group==1), Y(T.group==1), 'CNV SD', 'Alpha power', 'Young', 0, xlim_min)
scatter_and_corr(X(T.group==2), Y(T.group==2), 'CNV SD', 'Alpha power', 'Older', 0, xlim_min)

%%
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(Res(1:size(ERP_avg_amp_stdev_young, 2), 1), Res(1:size(ERP_avg_amp_stdev_young, 2), 2), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 1), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'ERP variability residuals')
% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(Res(1:size(ERP_avg_amp_stdev_young, 2), 2),  Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'CNV variability residuals')

%% regressout beta power from ERP variability of both groups
Res = [];
task = 1;
Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
X = [T.beta_power_simpleRT, ones(length(T.beta_power_simpleRT), 1)];
[B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
% scatter plot and correlation analyses
xlim_min = -.02;
scatter_and_corr(X(T.group==1), Y(T.group==1), 'CNV SD', 'Beta power', 'Young - simple RT', 0, xlim_min)
scatter_and_corr(X(T.group==2), Y(T.group==2), 'CNV SD', 'Beta power', 'Older - simple RT', 0, xlim_min)
task = 2;
Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
X = [T.beta_power_gng, ones(length(T.beta_power_simpleRT), 1)];
[B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
% scatter plot and correlation analyses - go/no-go
scatter_and_corr(X(T.group==1), Y(T.group==1), 'CNV SD', 'Beta power', 'Young', 0, xlim_min)
scatter_and_corr(X(T.group==2), Y(T.group==2), 'CNV SD', 'Beta power', 'Older', 0, xlim_min)


% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(Res(1:size(ERP_avg_amp_stdev_young, 2), 1), Res(1:size(ERP_avg_amp_stdev_young, 2), 2), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 1), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'ERP variability residuals')


%% regressout power of slow fluctuations as measuerd in raw PSD
% load PSD fooof variables for regression - FCz fitting with 4 peaks
% calculated in foof_analysis_power_spectrum_thresh1_4pks.m
table_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\';
load([table_dir, 'Frequencies.mat']); % load frequencies used in fooof analysis
load([table_dir, 'PowerSpectralDensity_Young']); load([table_dir, 'PowerSpectralDensity_Older']);
% avg power spectrum across trials - FCz - channel = 16
% group: 1 = young; 2 = older; task: 1  = simple RT; 2 = gng
chan = 16; clear avg_power_spectrum_young avg_power_spectrum_older
for p = 1:size(PowerSpectralDensity_Young, 1)
    avg_power_spectrum_young(1, p, :) = mean([PowerSpectralDensity_Young{p, 2, chan} PowerSpectralDensity_Young{p, 3, chan}], 2);
    avg_power_spectrum_young(2, p, :) = mean([PowerSpectralDensity_Young{p, 4, chan} PowerSpectralDensity_Young{p, 5, chan}], 2);
end

for p = 1:size(PowerSpectralDensity_Older, 1)
    avg_power_spectrum_older(1, p, :) =  mean([PowerSpectralDensity_Older{p, 2, chan} PowerSpectralDensity_Older{p, 3, chan}], 2);
    avg_power_spectrum_older(2, p, :) = mean([PowerSpectralDensity_Older{p, 4, chan} PowerSpectralDensity_Older{p, 5, chan}], 2);
end

% remove outliers
avg_power_spectrum_young(:, outliers_young, :) = [];
avg_power_spectrum_older(:, outliers_older, :) = [];

Res = [];
for task = 1:2

    Y = [squeeze(ERP_avg_amp_stdev_young(task, :, 16))'; squeeze(ERP_avg_amp_stdev_older(task, :, 16))'];
    X = [[squeeze(10*log10(avg_power_spectrum_young(task, :, 2)))'; squeeze(10*log10(avg_power_spectrum_older(task, :, 2)))'],...
        ones(size(avg_power_spectrum_young, 2)+size(avg_power_spectrum_older, 2), 1)];
    [B,BINT,Res(:, task),RINT,STATS] = regress(Y,X);
    if task == 2
        % scatter plot and correlation analyses
        %scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt, xlim_min)
        xlim_min = 9;
        scatter_and_corr(X(T.group==1), Y(T.group==1), 'CNV SD', 'log10(slow spectral power)', 'Young', 1, xlim_min)
        scatter_and_corr(X(T.group==2), Y(T.group==2), 'CNV SD', 'log10(low spectral power)', 'Older', 1, xlim_min)
    end

end

% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(Res(1:size(ERP_avg_amp_stdev_young, 2), 1), Res(1:size(ERP_avg_amp_stdev_young, 2), 2), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 1), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'CNV variability residuals')
% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(Res(1:size(ERP_avg_amp_stdev_young, 2), 2), Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2), 'CNV SD residuals')

[H,P,CI,STATS] = ttest2(Res(1:size(ERP_avg_amp_stdev_young, 2), 2),  Res(size(ERP_avg_amp_stdev_young, 2)+1:end, 2))
% T_values = [T_values, Res];
% 
% R = array2table(T_values, ...
%     'VariableNames',variable_names);
% 
% filename = 'ERPstdev_adjusted4exp_offset_slowpower.xlsx';
% writetable(R,filename,'Sheet',1,'Range','A1')

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
    axis([0 6 -7 10]);
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
    ax = gca; ax.FontSize = 24; ax.FontName = 'Arial';ax.Color = 'none';
    ylabel(ylabel_txt, 'FontSize', 32, 'FontWeight','normal');
    xlabel(xlabel_txt, 'FontSize', 32, 'FontWeight','normal');
    title(title_txt, 'FontSize', 32, 'FontWeight','normal');
    axis([xlim_min inf -inf inf])
    xl = xlim; yl = ylim;
    xt = .1*(xl(2)-xl(1))/4+xl(1); yt = 3.5*(yl(2)-yl(1))/4+yl(1);
    str = {['r = ' num2str(r(1, 2))], ['p = ' num2str(p(1, 2))]}
    if incl_txt
        text(xt,yt,str,'FontSize',14)
    end
end



function plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)

    % plot data for young group - go/nogo task
    figure; box on; hold on
    
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
    axis([0 3 -7 10]);
    ax = gca;
    c = ax.Color;
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