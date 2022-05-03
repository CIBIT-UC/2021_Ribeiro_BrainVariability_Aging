% correlations between spectral parameters from fooof analysis and CNV
% variability for all EEG channels
clear; close all
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% load spectral parameters
% calculated in C:\Users\Maria Ribeiro\Documents\GitHub\2021_Ribeiro_BrainVariability_Aging\Figure5\foof_analysis_power_spectrum_thresh1_4pks_all_chans.m
cd 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
load r_squared_error_older % (task, p, chan, r2 and error)
load r_squared_error_young
load exponent_young; load exponent_older; load offset_young; load offset_older; 
load oscillations_power_young; load oscillations_power_older

% load CNV variability data
% load ERP variability variables calculated in
% C:\Users\Maria Ribeiro\Documents\GitHub\2021_Ribeiro_BrainVariability_Aging\Figure2\Figure2C
load_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ERP_variability\';
load([load_dir 'ERP_std']); % group x task - participant x channel

% % find outliers - did not exclude outliers!
% for task = 1:3
%     for chan = 1:59
%         R2_transf_young(:, chan) = squeeze(fisherz(sqrt(r_squared_error_young(task, :, chan, 1))));
%         R2_transf_older(:, chan) = squeeze(fisherz(sqrt(r_squared_error_older(task, :, chan, 1))));      
%     end
%     
%     zscoreR2_young(task, :, :) = zscore(R2_transf_young, 0, 'all');
%     zscoreR2_older(task, :, :) = zscore(R2_transf_older, 0, 'all');
%     [index_part_older{task}, index_chan_older{task}] = find(zscore(R2_transf_older) < -3);
%     outliers_older{task} = unique(index_part_older{task});
%     [index_part_young{task}, index_chan_young{task}] = find(zscore(R2_transf_young) < -3);
%     outliers_young{task} = unique(index_part_young{task});
%     
% end
% % 
% % delete outliers
for task = 1:3
    exponent{1, task} = squeeze(exponent_young(task, :, :));
    exponent{2, task} = squeeze(exponent_older(task, :, :));
    
    offset{1, task} = squeeze(offset_young(task, :, :));
    offset{2, task} = squeeze(offset_older(task, :, :));
    
    oscillations_power{1, task} = squeeze(oscillations_power_young(task, :, :, :));
    oscillations_power{2, task} = squeeze(oscillations_power_older(task, :, :, :));
end
% 
% for task = 1:2 % only somple RT and gng
%     ERP_std{1, task}(outliers_young{task+1}, :) = [];% group x task - participant x channel
%     ERP_std{2, task}(outliers_older{task+1}, :) = [];% group x task - participant x channel
% end
% 
% for task = 1:3
%     exponent{1, task}(outliers_young{task}, :) = [];
%     exponent{2, task}(outliers_older{task}, :) = [];
%     
%     offset{1, task}(outliers_young{task}, :) = [];
%     offset{2, task}(outliers_older{task}, :) = [];
%     
%     oscillations_power{1, task}(outliers_young{task}, :, :) = []; % task x participants x channel x freq band (theta, alpha, beta)
%     oscillations_power{2, task}(outliers_older{task}, :, :) = [];
% end


% correlate spectral parameters with CNV variability
% PSD exponent
for grp = 1:2
    for task = 1:2 % simple RT, gng
%                 [pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX,dataY,n_perm,tail,alpha_level,stat,reports,seed_state)
% Required Inputs:
%  dataX   - A vector or 2D matrix of data (Observation x Variable)
%  dataY   - A vector or 2D matrix of data (Observation x Variable).  Each
%            element of dataX should be paired with an element of dataY.
        [pval(grp, task, :), corr_obs(grp, task, :), crit_corr(grp, task, :), est_alpha, seed_state] = mult_comp_perm_corr(ERP_std{grp, task},squeeze(exponent{grp, task+1}(:, :)));
    end
end

% plot corr coefficients
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
% title_txt = {'Simple RT' 'Go/no-go'};
for grp = 1:2
    for task = 2%1:2

        sig_r_values = zeros(1, 59);
        sig_r_values(squeeze(corr_obs(grp, task, :)) > crit_corr(grp, task, 2)) = 1;
        sig_r_values(squeeze(corr_obs(grp, task, :)) < crit_corr(grp, task, 1)) = 1;
        sig_chan_number = find(sig_r_values == 1);

        for x = 1:length(sig_chan_number)
            sig_chans_corr{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
        end


        figure;
        if ~isempty(sig_chan_number)
            topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
            'plotdisk', 'on') ; hold on % ,  'hcolor'  , 'none'
        end
        topoplot(corr_obs(grp, task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        Axis.FontSize = 16; caxis([-.65 .65]); 
        colorbar;
        colorbar('Ticks',[-.6 -.3 0 .3 .6], 'FontSize', 30, 'FontWeight','normal');
        colormap(crameri('tofino'));
    %     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
    %     set(get(gca,'title'),'Position',[0,-.65, 0])
    end
end


%% regress out the effect of exponent from CNV variability and compare across groups
% ERP_std = group x task - participant x channel
% exponent = group x task - participant x channel

for task = 1:2
    for chan = 1:59
        Y = [squeeze(ERP_std{1, task}(:, chan)); squeeze(ERP_std{2, task}(:, chan))];
        X = [[squeeze(exponent{1, task+1}(:, chan)); squeeze(exponent{2, task+1}(:, chan))] ones(length(Y), 1)];
        [~,~,Res(task, :, chan),~,~] = regress(Y,X);
    end
end

% compare residuals across groups
clear pval t_orig crit_t
for task = 1:2
    data1 = squeeze(Res(task, 1:size(ERP_std{1, task}, 1), :));
    data2 = squeeze(Res(task, size(ERP_std{1, task}, 1)+1:end, :));
    [pval(task, :), t_orig(task, :), crit_t(task, :), ~, ~]=mult_comp_perm_t2(data1,data2);
end

% plot t-values in scalp topography
clear sig_chans_corr
for task = 2%1:2
    sig_r_values = zeros(1, 59);
    sig_r_values(squeeze(t_orig(task, :)) > crit_t(task, 2)) = 1;
    sig_r_values(squeeze(t_orig(task, :)) < crit_t(task, 1)) = 1;
    sig_chan_number = find(sig_r_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans_corr{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end

    figure;
    if ~isempty(sig_chan_number)
        topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
        'plotdisk', 'on') ; hold on % ,  'hcolor'  , 'none'
    end
    topoplot(t_orig(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    Axis.FontSize = 16; caxis([0 4]); 
    colorbar;
    colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
    colormap(crameri('imola'));
    %     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
    %     set(get(gca,'title'),'Position',[0,-.65, 0])
end

%% correlation with PSD offset

% correlate spectral parameters with CNV variability
clear pval corr_obs crit_corr
for grp = 1:2
    for task = 1:2 % simple RT, gng
%                 [pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX,dataY,n_perm,tail,alpha_level,stat,reports,seed_state)
% Required Inputs:
%  dataX   - A vector or 2D matrix of data (Observation x Variable)
%  dataY   - A vector or 2D matrix of data (Observation x Variable).  Each
%            element of dataX should be paired with an element of dataY.
       [pval(grp, task, :), corr_obs(grp, task, :), crit_corr(grp, task, :), est_alpha, seed_state] = mult_comp_perm_corr(ERP_std{grp, task},squeeze(offset{grp, task+1}(:, :)));
    end
end

%% plot corr coefficients
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_txt = {'Simple RT' 'Go/no-go'};
for grp = 1:2
    for task = 2

        sig_r_values = zeros(1, 59);
        sig_r_values(squeeze(corr_obs(grp, task, :)) > crit_corr(grp, task, 2)) = 1;
        sig_r_values(squeeze(corr_obs(grp, task, :)) < crit_corr(grp, task, 1)) = 1;
        sig_chan_number = find(sig_r_values == 1);

        for x = 1:length(sig_chan_number)
            sig_chans_corr{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
        end


        figure;
        if ~isempty(sig_chan_number)
            topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
            'plotdisk', 'on') ; hold on % ,  'hcolor'  , 'none'
        end
        topoplot(corr_obs(grp, task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        Axis.FontSize = 16; caxis([-.8 .8]);
        colorbar;
        colorbar('Ticks',[-.8 -.4 0, .4, .8], 'FontSize', 30, 'FontWeight','normal');
        colormap(crameri('tofino'));
    %     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
    %     set(get(gca,'title'),'Position',[0,-.65, 0])
    end
end

% regress out the effect of offset from CNV variability and compare across groups
% ERP_std = group x task - participant x channel
% offset = group x task - participant x channel
clear Res
for task = 1:2
    for chan = 1:59
        Y = [squeeze(ERP_std{1, task}(:, chan)); squeeze(ERP_std{2, task}(:, chan))];
        X = [[squeeze(offset{1, task+1}(:, chan)); squeeze(offset{2, task+1}(:, chan))] ones(length(Y), 1)];
        [~,~,Res(task, :, chan),~,~] = regress(Y,X);
    end
end

% compare residuals across groups
clear pval t_orig crit_t
for task = 1:2
    data1 = squeeze(Res(task, 1:size(ERP_std{1, task}, 1), :));
    data2 = squeeze(Res(task, size(ERP_std{1, task}, 1)+1:end, :));
    [pval(task, :), t_orig(task, :), crit_t(task, :), ~, ~]=mult_comp_perm_t2(data1,data2);
end

% plot t-values in scalp topography
for task = 2%1:2
    sig_r_values = zeros(1, 59);
    sig_r_values(squeeze(t_orig(task, :)) > crit_t(task, 2)) = 1;
    sig_r_values(squeeze(t_orig(task, :)) < crit_t(task, 1)) = 1;
    sig_chan_number = find(sig_r_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans_corr{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end

    figure;
    if ~isempty(sig_chan_number)
        topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
        'plotdisk', 'on') ; hold on % ,  'hcolor'  , 'none'
    end
    topoplot(t_orig(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    Axis.FontSize = 16; caxis([0 4]); 
    colorbar;
    colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
    colormap(crameri('imola'));
    %     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
    %     set(get(gca,'title'),'Position',[0,-.65, 0])
end


%% slow PSD power correlation with ERP variability
% calculated in \2021_Ribeiro_BrainVariability_Aging\Figure5\pre_stimulus_power_spectrum.m
table_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\three_seconds_baseline_power_spectra\';
load([table_dir, 'Frequencies.mat']); % load frequencies used in fooof analysis
load([table_dir, 'PowerSpectralDensity_Young']); load([table_dir, 'PowerSpectralDensity_Older']); %{participant x run x channel}
grp = 2;
for p = 1:size(PowerSpectralDensity_Older, 1)
    for chan = 1:size(PowerSpectralDensity_Older, 3)
        % gng task = runs 3 and 4
%         figure; plot(Frequencies, mean([PowerSpectralDensity_Older{1, 3, chan}, PowerSpectralDensity_Older{1, 4, chan}], 2))
        slow_fluctuations{grp}(p, chan) = log10(mean([PowerSpectralDensity_Older{p, 3, chan}(2, :), PowerSpectralDensity_Older{p, 4, chan}(2, :)], 2));
    end
end

grp = 1;
for p = 1:size(PowerSpectralDensity_Young, 1)
    for chan = 1:size(PowerSpectralDensity_Young, 3)
        % gng task = runs 3 and 4
%         figure; plot(Frequencies, mean([PowerSpectralDensity_Older{1, 3, chan}, PowerSpectralDensity_Older{1, 4, chan}], 2))
        slow_fluctuations{grp}(p, chan) = log10(mean([PowerSpectralDensity_Young{p, 3, chan}(2, :), PowerSpectralDensity_Young{p, 4, chan}(2, :)], 2));
    end
end

% correlate amplitude of slow fluctuations with CNV variability
% plot scatter plot for channel FCz
for ch = 1:length(chanlocs_EEGChanOnly)
    if strcmp(chanlocs_EEGChanOnly(ch).labels, 'FCz')
        fcz_chan = ch;
    end
end
clear pval corr_obs crit_corr
for grp = 1:2
    for task = 2 % simple RT, gng
%                 [pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX,dataY,n_perm,tail,alpha_level,stat,reports,seed_state)
% Required Inputs:
%  dataX   - A vector or 2D matrix of data (Observation x Variable)
%  dataY   - A vector or 2D matrix of data (Observation x Variable).  Each
%            element of dataX should be paired with an element of dataY.
        [pval(grp, task, :), corr_obs(grp, task, :), crit_corr(grp, task, :), est_alpha, seed_state] = mult_comp_perm_corr(ERP_std{grp, task},slow_fluctuations{grp}(:, :));
%         scatter_and_corr(X, Y, ylabel_txt, xlabel_txt, title_txt, incl_txt, xlim_min)
%         scatter_and_corr(slow_fluctuations{grp}(:, fcz_chan), ERP_std{grp, task}(:, fcz_chan), ylabel_txt, xlabel_txt, title_txt, incl_txt, xlim_min)
    end
end

% plot corr coefficients
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
% title_txt = {'Simple RT' 'Go/no-go'};
for grp = 1:2
    for task = 2%1:2

        sig_r_values = zeros(1, 59);
        sig_r_values(squeeze(corr_obs(grp, task, :)) > crit_corr(grp, task, 2)) = 1;
        sig_r_values(squeeze(corr_obs(grp, task, :)) < crit_corr(grp, task, 1)) = 1;
        sig_chan_number = find(sig_r_values == 1);

        for x = 1:length(sig_chan_number)
            sig_chans_corr{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
        end


        figure;
        if ~isempty(sig_chan_number)
            topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
            'plotdisk', 'on') ; hold on % ,  'hcolor'  , 'none'
        end
        topoplot(corr_obs(grp, task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        Axis.FontSize = 16; caxis([-1 1]); 
        colorbar;
        colorbar('Ticks',[-1 -.5 0 .5 1], 'FontSize', 30, 'FontWeight','normal');
        colormap(crameri('tofino'));
    %     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
    %     set(get(gca,'title'),'Position',[0,-.65, 0])
    end
end

% regress out the effect of slow signal fluctuations from CNV variability and compare across groups
% ERP_std = group x task - participant x channel
% offset = group x task - participant x channel
clear Res
for chan = 1:59
    Y = [squeeze(ERP_std{1, task}(:, chan)); squeeze(ERP_std{2, task}(:, chan))];
    X = [[squeeze(slow_fluctuations{1}(:, chan)); squeeze(slow_fluctuations{2}(:, chan))] ones(length(Y), 1)];
    [~,~,Res(:, chan),~,~] = regress(Y,X);
end

% compare residuals across groups
clear pval t_orig crit_t
data1 = squeeze(Res(1:size(ERP_std{1, task}, 1), :));
data2 = squeeze(Res(size(ERP_std{1, task}, 1)+1:end, :));
[pval, t_orig, crit_t, ~, ~]=mult_comp_perm_t2(data1,data2);

% plot t-values in scalp topography
sig_r_values = zeros(1, 59);
sig_r_values(squeeze(t_orig) > crit_t(2)) = 1;
sig_r_values(squeeze(t_orig) < crit_t(1)) = 1;
sig_chan_number = find(sig_r_values == 1);

for x = 1:length(sig_chan_number)
    sig_chans_corr{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
end

figure;
if ~isempty(sig_chan_number)
    topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on') ; hold on % ,  'hcolor'  , 'none'
end
topoplot(t_orig, chanlocs_EEGChanOnly, 'electrodes', 'off'); 
Axis.FontSize = 16; caxis([0 4]); 
colorbar;
colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
colormap(crameri('imola'));
%     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
%     set(get(gca,'title'),'Position',[0,-.65, 0])

%% correlation between alpha and beta oscillations and ERP variability
% calculated in C:\Users\Maria Ribeiro\Documents\GitHub\2021_Ribeiro_BrainVariability_Aging\Figure5\foof_analysis_power_spectrum_thresh1_4pks_all_chans.m
% oscillations_power_older(task, p, chan, 1) = task x participants x channel x freq band (theta, alpha, beta)

close all
% correlate alpha and beta power with CNV variability
for osc = 2:3 %alpha and beta peaks
    clear pval corr_obs crit_corr
    for grp = 1:2
        for task = 1:2 % simple RT, gng
        % [pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX,dataY,n_perm,tail,alpha_level,stat,reports,seed_state)
        % Required Inputs:
        % dataX   - A vector or 2D matrix of data (Observation x Variable)
        % dataY   - A vector or 2D matrix of data (Observation x Variable).  Each element of dataX should be paired with an element of dataY.
           [pval(grp, task, :), corr_obs(grp, task, :), crit_corr(grp, task, :), est_alpha, seed_state] = mult_comp_perm_corr(ERP_std{grp, task},squeeze(oscillations_power{grp, task+1}(:, :, osc)));
        end
    end

    % plot corr coefficients
    load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
    title_txt = {'Simple RT' 'Go/no-go'};
    for grp = 1:2
        for task = 2

            sig_r_values = zeros(1, 59);
            sig_r_values(squeeze(corr_obs(grp, task, :)) > crit_corr(grp, task, 2)) = 1;
            sig_r_values(squeeze(corr_obs(grp, task, :)) < crit_corr(grp, task, 1)) = 1;
            sig_chan_number = find(sig_r_values == 1);

            for x = 1:length(sig_chan_number)
                sig_chans_corr{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
            end


            figure;
            if ~isempty(sig_chan_number)
                topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                'plotdisk', 'on') ; hold on % ,  'hcolor'  , 'none'
            end
            topoplot(corr_obs(grp, task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
            if osc == 2
                Axis.FontSize = 16; caxis([-.6 .6]);
            else
                Axis.FontSize = 16; caxis([-.5 .5]);
            end
            colorbar;
            colorbar('Ticks',[-.4 0, .4], 'FontSize', 30, 'FontWeight','normal');
            colormap(crameri('tofino'));
        %     title(title_txt{task}, 'FontSize', 24, 'FontWeight','normal')
        %     set(get(gca,'title'),'Position',[0,-.65, 0])
        
        
        end
    end
end

%% compare CNV variability across groups after adjusting for slow signal fluctuations





%% 
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