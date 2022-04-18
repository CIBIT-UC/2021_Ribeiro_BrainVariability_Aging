% relies on 
% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
% https://www.mathworks.com/matlabcentral/fileexchange/54585-mult_comp_perm_t2-data1-data2-n_perm-tail-alpha_level-mu-t_stat-reports-seed_state?s_tid=srchtitle

% change eeglab function topoplot so that electrode disks are plotted in
% black - line 275 COLORARRAY  = { [0 0 0] [0.5 0 0] [0 0 0] };

clear; close all;
% participants id
% eeg included participants
young_eeg=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older_eeg=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];

group = {young_eeg, older_eeg};
task={'D1', 'D2', 'G1', 'G2'};
number_of_trials_eeg = {};
% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

betas_young = cell(2, 1); betas_older = cell(2, 1); r_young = cell(2, 1); r_older = cell(2, 1);
for grp = 1:2
    part = 0;
    for p = group{grp}
        part = part + 1;
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
             
        % erp data - dependent variable
        eeg_directory = strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\EEG\');

        % create variable with data electrodes X time frames X trials -
        % starting at cue onset
        ERPs= {}; RT_incl = {}; number_of_trials = [];
       for t = 1:4
           % clear eeglab
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            % load file calculated in G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\eeglab_analysis_3.m
           filename=strcat('AB', num2str(p), '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej_epochs_baselined.set');
           EEG = pop_loadset(filename, eeg_directory);
           
           % interpolate bad channels removed earlier
            if length(EEG.chanlocs)<67 && p ~= 6 && p~=38 % these participants do not have eyetracking data but have all channels ok
                [EEG] = interpol( EEG);
            end
           
                  
           % calculate reaction time for every trial and error trials to
           % exclude
            clear epoch_event;
            epoch_event=zeros(size(EEG.event, 2), 3);

            for e=1:size(EEG.event, 2)
                epoch_event(e, 1)=EEG.event(e).epoch;
                epoch_event(e, 2)=str2double(EEG.event(e).type);
                epoch_event(e, 3)=EEG.event(e).latency;
            end

            cue_only_correct_trial=[];
            correct_trial=[];
            response2cue=[];
            multiple_responses=[];
            misses=[];
            response2cue=[];
            RT=[];
            slow_responses=[];
            slow_RT=[];
            cue2target_interval = [];
            response2nogo=[];

            for r = 1:epoch_event(end, 1)
                events_per_epoch_index=[];
                events_per_epoch_index=find(epoch_event(:,1)==r);

                events_per_epoch=epoch_event(events_per_epoch_index, 2);
                latencies_per_epoch=epoch_event(events_per_epoch_index, 3);

                if length(events_per_epoch)==1
                    cue_only_correct_trial=[cue_only_correct_trial; r];
                else
                    if events_per_epoch(2)==2 && events_per_epoch(3)==5 && length(events_per_epoch)==3
                        correct_trial=[correct_trial; r];
                        RT=[RT; r (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/500]; 
                        cue2target_interval = [cue2target_interval; r (latencies_per_epoch(2)-latencies_per_epoch(1))*1000/500];
                    elseif events_per_epoch(2)==2 && events_per_epoch(3)==5 && events_per_epoch(4)==3 
                        slow_responses=[slow_responses; r];
                        slow_RT=[slow_RT; r (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/500];
                    elseif events_per_epoch(2)==5
                        response2cue=[response2cue; r];
                    elseif length(events_per_epoch)>3 && events_per_epoch(2)==2 && events_per_epoch(3)==5 && events_per_epoch(4)==5
                        multiple_responses=[multiple_responses; r]; 
                    elseif events_per_epoch(2)==2 && events_per_epoch(3)==3
                        misses=[misses; r];
                    elseif length(events_per_epoch)>=3 && (events_per_epoch(2)==4 && events_per_epoch(3)==6 && t>2)
                        response2nogo=[response2nogo; r];  
                    end
                end 
            end

            % exclude error trials and trials after error - there is a problem
            % where chunks of EEG with muscle artifact were removed trials
            % after error might NOT immediately follow an error - shoudl check
            % with urevent
            % error_trials - all errors, response2cue, multiple responses, misses, slow responses
            error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses, response2nogo);
            % only correct go trials, excluding trials after error
            include_epochs = setdiff(correct_trial, error_trials+1);
            number_of_trials(t) = length(include_epochs);
            EEG = pop_select( EEG, 'trial',include_epochs);
            cnt = 0;
            for inc = include_epochs'
                cnt = cnt + 1;
                RT_incl{t}(cnt, 1) = RT(RT(:, 1) == inc, 2);
    %             cue2target_interval_inc{t}(cnt, 1) = cue2target_interval(cue2target_interval(:, 1) == inc, 2);
            end
        
           % ERP data - avg amplitude between 1000 and 1500 ms after cue
           % onset
           ERPs{t} = squeeze(mean(EEG.data(1:59, 600:850, :), 2)); % chan X trials       
       end 
       
       number_of_trials_eeg{grp}(part, :) = [number_of_trials(1)+number_of_trials(2), number_of_trials(3)+number_of_trials(4)];
         %% create variables ERP ampltiude and ERP std for each channel

         % tasks: simple RT (runs 1 and 2) and gng (runs 3 and 4)
         ERP_amp{grp, 1}(part, :) = mean([ERPs{1}, ERPs{2}], 2);
         ERP_amp{grp, 2}(part, :) = mean([ERPs{3}, ERPs{4}], 2);
         ERP_std{grp, 1}(part, :) = std([ERPs{1}, ERPs{2}],[], 2);
         ERP_std{grp, 2}(part, :) = std([ERPs{3}, ERPs{4}],[], 2);
    
    end
end
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability');
save number_of_trials_eeg number_of_trials_eeg
save ERP_amp ERP_amp
save ERP_std ERP_std

% compare number of trials across groups
for task =1:2
    [h,p,ci,stats] = ttest2(number_of_trials_eeg{1}(:, task), number_of_trials_eeg{2}(:, task))
    mean_trials(1, task) = mean(number_of_trials_eeg{1}(:, task));
    sd_trials(1, task) = std(number_of_trials_eeg{1}(:, task));
    mean_trials(2, task) = mean(number_of_trials_eeg{2}(:, task));
    sd_trials(2, task) = std(number_of_trials_eeg{2}(:, task));
end
%% compare CNV amplitude and variability across groups
% stats using permutation method
% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability');
load ERP_amp; load ERP_std
for task = 1:2
    [pval_amp(task, :), t_orig_amp(task, :), crit_t_amp(task, :), est_alpha, seed_state]=mult_comp_perm_t2(ERP_amp{1, task},ERP_amp{2, task});
    [pval_std(task, :), t_orig_std(task, :), crit_t_std(task, :), est_alpha, seed_state]=mult_comp_perm_t2(ERP_std{1, task},ERP_std{2, task});
end


% plot  - t values
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_txt = {'Simple RT' 'Go/no-go'};
for task = 1:2
    
    sig_t_values = zeros(1, 59);
    sig_t_values(t_orig_std(task, :) > crit_t_std(task, 2)) = 1;
    sig_t_values(t_orig_std(task, :) < crit_t_std(task, 1)) = 1;
    sig_chan_number = find(sig_t_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans_std{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end
    
    
    figure;
    topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on') ; hold on % ,  'hcolor'  , 'none'
    topoplot(t_orig_std(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    Axis.FontSize = 16; caxis([0 6]); 
    colorbar;
    colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
    colormap(crameri('imola'));
%     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
%     set(get(gca,'title'),'Position',[0,-.65, 0])

end

%%
for task = 1:2
    
    sig_t_values = zeros(1, 59);
    sig_t_values(t_orig_amp(task, :) > crit_t_amp(task, 2)) = 1;
    sig_t_values(t_orig_amp(task, :) < crit_t_amp(task, 1)) = 1;
    sig_chan_number = find(sig_t_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans_amp{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end
    
    figure; 
    topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
        'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
    topoplot(t_orig_amp(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off');
    Axis.FontSize = 16; caxis([0 5]); 
    colorbar;
    colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
    colormap(crameri('imola'));
%     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
%     set(get(gca,'title'),'Position',[0,-.65, 0])
%     hold on


end


%% plot
clear sig_chans
load ERP_amp; load ERP_std
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat;
title_txt = {'Young - simple RT' 'Older - simple RT'; 'Young - go/no-go' 'Older - go/no-go'};
for grp = 1:2
    for task = 1:2
        figure;
        topoplot(mean(ERP_amp{grp, task}, 1), chanlocs_EEGChanOnly); 
        caxis([-10 0]); c.Axis.FontSize = 16;
        colorbar;
        colorbar('Ticks',[-10, -5, 0], 'FontSize', 30, 'FontWeight','normal');
        colormap(crameri('imola'));
%         title(title_txt{task, grp}, 'FontSize', 24, 'FontWeight','normal')
%         set(get(gca,'title'),'Position',[0,-1, 0])
        
        figure;
        topoplot(mean(ERP_std{grp, task}, 1), chanlocs_EEGChanOnly); 
        caxis([0 12]); c.Axis.FontSize = 16;
        colorbar;
        colorbar('Ticks',[0, 5, 10], 'FontSize',30, 'FontWeight','normal');
        colormap(crameri('imola'));
%         title(title_txt{task, grp}, 'FontSize', 24, 'FontWeight','normal')
%         set(get(gca,'title'),'Position',[0,-1, 0])
        
    end
end

%%
% plot_data_all_electrodes(coeff_robust_chan{1, 1}, coeff_robust_chan{2, 1}, 'Amplitude', 'Simple RT')
% 
% plot_data_all_electrodes(coeff_robust_chan{1, 2}, coeff_robust_chan{2, 2}, 'Amplitude', 'Go/no-go')

%% graphs of coefficients
% % FCz CNV amplitude; phasic pupil; cue-to-target interval; time-on-task - total 4 continuous variables
%     % interactions ERP x task, pupil x task, time-on-task x task
% y_label_text = 'Estimated coefficients';
% task_name = {'Simple RT' ' Go/no-go'};
% title_text = {'CNV amplitude' 'Pupil response' 'Cue-to-target interval'}; 
% 
% for task = 1:2 
%     for coef = 1:3
%         % plot_all_data_2groups(data_grp1, data_grp2, y_label_text)
%         plot_all_data_2groups(betas_young{task}(: , coef+1), betas_older{task}(: , coef+1), y_label_text, [title_text(coef), task_name(task)])
%     end
% end

%% function plot data 2 groups 1 task
function plot_all_data_2groups(data_grp1, data_grp2, y_label_text, title_text)
       figure;
    % plot data for group 1
        yMean=nanmean(data_grp1);
    y_se = std(data_grp1)/sqrt(length(data_grp1));
    
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
    %plot line at zero    
    plot([0 4],[0 0],'--k','LineWidth', 1);
     %plot the mean line    
    plot([1-0.3 1+0.3],[yMean yMean] ,'Color','k','LineWidth',5);
   

    %% group 2
    yMean=nanmean(data_grp2);
    y_se = std(data_grp2)/sqrt(length(data_grp2));
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
    ax.FontSize = 22;
    ax.FontName = 'Arial';
    ax.Color = 'none';
%     ax.XTickLabel= [];
    ax.XTickLabel= [{'' 'Young' 'Older' ''}];
%     xticks([1 2])
    ylabel(y_label_text, 'FontSize', 24, 'FontWeight','normal')
    title(title_text, 'FontSize', 24, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end



function plot_data_all_electrodes(data_young, data_older, yaxis_title, title_txt)

    % plot 
    figure; 
    plot(.8:58.8, data_young, 'o', 'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k', 'MarkerSize', 5); hold on;
    plot(1.2:59.2, data_older, 'o', 'MarkerFaceColor',[1 .8 .8], 'MarkerEdgeColor','r','MarkerSize', 5); hold on
    mean_data = mean(data_young, 1, 'omitnan'); se_data = std(data_young, [], 1, 'omitnan')/sqrt(size(data_young, 1));
    plot(1:59, mean_data', '-', 'color', 'k', 'LineWidth', 2); hold on;
    jbfill(1:59, mean_data+se_data, mean_data-se_data,'k','k', 0.1); hold on;

    mean_data = mean(data_older, 1, 'omitnan'); se_data = std(data_older, [], 1, 'omitnan')/sqrt(size(data_older, 1));
    plot(1:59,  mean_data', '-', 'color', 'r', 'LineWidth', 2, 'MarkerEdgeColor', 'r',...
        'MarkerSize', 10); hold on;
    jbfill(1:59, mean_data+se_data, mean_data-se_data,'r','r', 0.1); hold on;
    
    % plot zero line
    plot(0:60, zeros(1,61), '--k'); hold on
    

    axis([0 61 -inf inf]);
    ax = gca; c = ax.Color; ax.FontSize = 12; ax.FontName = 'Arial'; ax.Color = 'none';
%     set(gca, 'YScale', 'log')
%     ax.XTickLabel=x_axis(1:59);
%     xticks(1:59)
    xlabel('Electrodes', 'FontSize', 18, 'FontWeight','normal')
    ylabel(yaxis_title, 'FontSize', 18, 'FontWeight','normal')
    title(title_txt, 'FontSize', 18, 'FontWeight','normal')
    x0=0; y0=0; width=1600; height=300; set(gcf,'position',[x0,y0,width,height])
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
