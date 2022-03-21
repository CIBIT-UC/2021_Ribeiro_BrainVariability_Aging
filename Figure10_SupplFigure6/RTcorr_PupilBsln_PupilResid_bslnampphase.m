% across trials standard deviation of amplitude of cue-locked pupil
% response durign the preparatory period before target onset
clear; close all;
younger=[4 9 10 13 15 16 25 26 28 31 33 34 36 42 44 45 46 50 51 53 54 56 59 62 66 68 72 74 76 78 80 81 82 84 85];
older=[7 8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77  79 83   86];
task={'D1', 'D2', 'G1', 'G2'};
group = {younger older};
% std_pupil = cell(2, 5); std_pupil_avg_response = cell(2, 2); pupil_avg_response = cell(2, 2);
% pupil_cuelocked_response = cell(2, 2);
% pupil_percent_signal_change = cell(2, 5); 
% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for grp = 1:2
    count = 0;
    for p = group{grp}
        count = count + 1;
        pupil_directory=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\');
        load([pupil_directory, 'RejectEpochsCorrect']); % data from G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\pupil_analysis_3b_RejectedEpochs_SaveRunData.m
        
%         if intersect(p, [8, 11, 14, 15]) % in these participants use second repetition of recordign due to technical problems
%             task = {'W2', 'D1', 'D2', 'G1', 'G2'};
%         else
%             task={'W1', 'D1', 'D2', 'G1', 'G2'};
%         end
        
        clear pupil_avg_response_tmp pupil_response
        RT = cell(4, 1); time_on_task = cell(4, 1);
        for t = 1:4 % task runs excluding first run = passive task
        
            % calculate pupil go trials data
            filename=strcat(pupil_directory, 'AB', num2str(p), '_', task(t), '_PupilData.mat');
            load(filename{1})
            % use smoothed pupil diameter data with blink artefacts interpolated
            % measured in horizontal direction (x)
            pupil_data_percent_signal_change = PupilData(:, 3)/mean(PupilData(:, 3)) - 1;
            
%             save([pupil_directory, 'pupil_data_percent_signal_change.mat'], 'pupil_data_percent_signal_change');
            pupil_data_percent_signal_change(:, 2) = PupilData(:, 7); % add triggers data
%             filename = [pupil_directory, 'pupil_data_percent_signal_change.mat'];
            setname = ['AB', num2str(p), '_', task{t}];
            EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',pupil_data_percent_signal_change','setname',setname,'srate',240);%,'pnts',0,'xmin',0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
            EEG = pop_chanevent(EEG, 2,'edge','leading','edgelen',0);
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = pop_epoch( EEG, {  '1'  }, [-.2  6]); % epoch locked to cue onset
            % delete epochs with artefacts as well as cue-only trials and error trials and trials after error trials
            EEG = pop_selectevent( EEG, 'omitepoch', RejectEpochsCorrect{t} ,'deleteevents','off','deleteepochs','on','invertepochs','off');
            pupil_baseline{t} = squeeze(mean(EEG.data(1, 1:48, :), 2));
            EEG = pop_rmbase( EEG, [-200    0]); % subtract baseline pupil
            
            % average pupil dilation from 1000 - 1500 ms after cue onset -
            % 240 sampling rate
            % pupil_avg_response{t} = squeeze(mean(EEG.data(1, 360:end, setdiff(1:60, RejectEpochsCorrect{t-1})), 2)); % time x trials 
            pupil_avg_response{t} = squeeze(mean(EEG.data(1, 481:600, :), 2)); % time x trials 
            
            % determine reaction time for each trial
            epoch_event=zeros(size(EEG.event, 2), 3);

            for e=1:size(EEG.event, 2)
                epoch_event(e, 1)=EEG.event(e).epoch;
                epoch_event(e, 2)=double(EEG.event(e).type);
                epoch_event(e, 3)=EEG.event(e).latency;
            end
            
            for r = 1:epoch_event(end, 1)
                events_per_epoch_index=[];
                events_per_epoch_index=find(epoch_event(:,1)==r);
                events_per_epoch=epoch_event(events_per_epoch_index, 2);
                latencies_per_epoch=epoch_event(events_per_epoch_index, 3);

                if events_per_epoch(2)==2 && events_per_epoch(3)==5 && length(events_per_epoch)==3
                    RT{t}=[RT{t}; (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/240]; 
%                     cue2target_interval = [cue2target_interval;(latencies_per_epoch(2)-latencies_per_epoch(1))*1000/240];
                end 
            end
            % create variable with trial number for each run - time on task
            time_on_task{t}(:,1)=setdiff(1:60, RejectEpochsCorrect{t})';% exclude trials that had artefacts

        end
        
                % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        % calculate phase of pupil fluctuations at cue-onset
        for t = 1:4 % task runs excluding first run = passive task
                % calculate pupil go trials data
                filename=strcat(pupil_directory, 'AB', num2str(p), '_', task(t), '_PupilData.mat');
                load(filename{1})
                % use smoothed pupil diameter data with blink artefacts interpolated
                % measured in horizontal direction (x)
                pupil_data_percent_signal_change = PupilData(:, 3)/mean(PupilData(:, 3));
                pupil_data_percent_signal_change(:, 2) = PupilData(:, 7); % add triggers data
                save([pupil_directory, 'pupil_data_percent_signal_change.mat'], 'pupil_data_percent_signal_change');
                filename = [pupil_directory, 'pupil_data_percent_signal_change.mat'];
                setname = ['AB', num2str(p), '_', task{t}];
                EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',filename,'setname',setname,'srate',240);%,'pnts',0,'xmin',0);
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
                EEG = pop_chanevent(EEG, 2,'edge','leading','edgelen',0);
                % filter data between .1 and .9 Hz
                EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'hicutoff',0.9);
                % hilbert transform
                y = hilbert(EEG.data);
                signal_phase = angle(y); % extract instantaneous phase
                signal_amplitude = abs(y); % extract amplitude envelop
                EEG.data(1, :) = signal_phase;
                EEG.data(2, :) = signal_amplitude;
                EEG = pop_epoch( EEG, {  '1'  }, [-.2  1.5]); % epoch locked to cue onset
                pupil_bsln_phase{t} = squeeze(EEG.data(1, 48, setdiff(1:60, RejectEpochsCorrect{t}))); % phase at cue-onset for each trial
                pupil_bsln_amp{t} = squeeze(EEG.data(2, 48, setdiff(1:60, RejectEpochsCorrect{t}))); % amplitude envelop at cue-onset for each trial
        end
        
        
%         % within-subject correlation between pupil response amplitude and reaction time 
%         [r,~,~,~,~,~] = skipped_correlation([pupil_avg_response{1}; pupil_avg_response{2}], [RT{1}; RT{2}],0);
%         coeff_robust_TPR_RT{grp, 1}(count) = r.Pearson;
%         
%         [r,~,~,~,~,~] = skipped_correlation([pupil_avg_response{3}; pupil_avg_response{4}], [RT{3}; RT{4}],0);
%         coeff_robust_TPR_RT{grp, 2}(count) = r.Pearson;
        
       % regress phase of slow fluctuations out of erp amplitude
       % simple RT
       y = [pupil_avg_response{1}; pupil_avg_response{2}];
       x = [[pupil_bsln_amp{1}.*cos(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*cos(pupil_bsln_phase{2})], [pupil_bsln_amp{1}.*sin(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*sin(pupil_bsln_phase{2})], ones(length(y), 1)];
       [~,~,Res_simpleRT] = regress(y, x);

       % gng
       y = [pupil_avg_response{3}; pupil_avg_response{4}];
       x = [[pupil_bsln_amp{3}.*cos(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*cos(pupil_bsln_phase{4})], [pupil_bsln_amp{3}.*sin(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*sin(pupil_bsln_phase{4})], ones(length(y), 1)];
       [~,~,Res_gng] = regress(y, x);
        
        % within-subject correlation between pupil response amplitude residuals and reaction time 
        [r,~,~,~,~,~] = skipped_correlation(Res_simpleRT, [RT{1}; RT{2}],0);
        coeff_robust_pupilresidphase_RT{grp, 1}(count) = r.Pearson;
        
        [r,~,~,~,~,~] = skipped_correlation(Res_gng, [RT{3}; RT{4}],0);
        coeff_robust_pupilresidphase_RT{grp, 2}(count) = r.Pearson;
        
        
        % within-subject correlation between pupil phase (cos and sin) and
        % RT
         y = [RT{1}; RT{2}];
         x = [[pupil_bsln_amp{1}.*cos(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*cos(pupil_bsln_phase{2})], [pupil_bsln_amp{1}.*sin(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*sin(pupil_bsln_phase{2})], ones(length(y), 1)];
         coefs_phase_RT{grp, 1}(count, :) = regress(y, x);
         y = [RT{3}; RT{4}];
         x = [[pupil_bsln_amp{3}.*cos(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*cos(pupil_bsln_phase{4})], [pupil_bsln_amp{3}.*sin(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*sin(pupil_bsln_phase{4})], ones(length(y), 1)];
         coefs_phase_RT{grp, 2}(count, :) = regress(y, x);
         
         
        % within-subject correlation between pupil phase (cos and sin) and
        % RT - controling for time-on-task - see 2016 van den brink et al Pupil diameter tracks lapses of attention
         y = [RT{1}; RT{2}];
         x = [[pupil_bsln_amp{1}.*cos(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*cos(pupil_bsln_phase{2})], [pupil_bsln_amp{1}.*sin(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*sin(pupil_bsln_phase{2})],...
             [time_on_task{1}; time_on_task{2}], ones(length(y), 1)];
         coefs_phase_timeontask_RT{grp, 1}(count, :) = regress(y, x);
         y = [RT{3}; RT{4}];
         x = [[pupil_bsln_amp{3}.*cos(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*cos(pupil_bsln_phase{4})], [pupil_bsln_amp{3}.*sin(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*sin(pupil_bsln_phase{4})],...
             [time_on_task{3}; time_on_task{4}], ones(length(y), 1)];
         coefs_phase_timeontask_RT{grp, 2}(count, :) = regress(y, x);
%         
%         
%          % within-subject correlation between pupil baseline and TPR amplitude
%         [r,~,~,~,~,~] = skipped_correlation([pupil_baseline{1}; pupil_baseline{2}], [pupil_avg_response{1}; pupil_avg_response{2}],0);
%         coeff_robust_pupilBsln_TPR{grp, 1}(count) = r.Pearson;
%         
%         [r,~,~,~,~,~] = skipped_correlation([pupil_baseline{3}; pupil_baseline{4}], [pupil_avg_response{3}; pupil_avg_response{4}],0);
%         coeff_robust_pupilBsln_TPR{grp, 2}(count) = r.Pearson;
        
        
        % within-subject correlation between pupil baseline and TPR
        % amplitude residuals
        [r,~,~,~,~,~] = skipped_correlation([pupil_baseline{1}; pupil_baseline{2}], Res_simpleRT, 0);
        coeff_robust_pupilBsln_TPRresid{grp, 1}(count) = r.Pearson;
        
        [r,~,~,~,~,~] = skipped_correlation([pupil_baseline{3}; pupil_baseline{4}], Res_gng, 0);
        coeff_robust_pupilBsln_TPRresid{grp, 2}(count) = r.Pearson;
        
        
%         % within-subject correlation between pupil baseline and RT
%         [r,~,~,~,~,~] = skipped_correlation([pupil_baseline{1}; pupil_baseline{2}], [RT{1}; RT{2}], 0);
%         coeff_robust_pupilBsln_RT{grp, 1}(count) = r.Pearson;
%         
%         [r,~,~,~,~,~] = skipped_correlation([pupil_baseline{3}; pupil_baseline{4}], [RT{3}; RT{4}], 0);
%         coeff_robust_pupilBsln_RT{grp, 2}(count) = r.Pearson;
        
        
%         % within-subject correlation between pupil baseline and cosine of
%         % baseline angle
%         [r,~,~,~,~,~] = skipped_correlation([pupil_baseline{1}; pupil_baseline{2}], [cos(pupil_bsln_phase{1}); cos(pupil_bsln_phase{2})], 0);
%         coeff_robust_pupilBsln_cos{grp, 1}(count) = r.Pearson;
%         
%         [r,~,~,~,~,~] = skipped_correlation([pupil_baseline{3}; pupil_baseline{4}], [cos(pupil_bsln_phase{3}); cos(pupil_bsln_phase{4})], 0);
%         coeff_robust_pupilBsln_cos{grp, 2}(count) = r.Pearson;

        
        % calculate std across trials, simple RT and gng task - pupil
        % residuals after reggressing pupil baseline
        pupil_resd_phase_std{grp, 1}(count) = std(Res_simpleRT); % simple RT
        pupil_resd_phase_std{grp, 2}(count) = std(Res_gng); % gng
        
        
        

    end
end

%% save variable
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability')
% save coeff_robust_pupilBsln_RT coeff_robust_pupilBsln_RT
% save coeff_robust_pupilBsln_cos coeff_robust_pupilBsln_cos
save coefs_phase_timeontask_RT coefs_phase_timeontask_RT
save pupil_resd_phase_std pupil_resd_phase_std
save coeff_robust_pupilBsln_TPRresid coeff_robust_pupilBsln_TPRresid
% % save coeff_robust_pupilBsln_TPR coeff_robust_pupilBsln_TPR
save coeff_robust_pupilresidphase_RT coeff_robust_pupilresidphase_RT % correlation coefficients RT vs pupil adjusted for basln phase
save coefs_phase_RT coefs_phase_RT
% save coeff_robust_TPR_RT coeff_robust_TPR_RT


%% within-subject correlation between pupil baseline and reaction time 
load coeff_robust_pupilBsln_RT
plot_all_data_2tasks(coeff_robust_pupilBsln_RT{1,1}, coeff_robust_pupilBsln_RT{1, 2},...
    coeff_robust_pupilBsln_RT{2, 1}, coeff_robust_pupilBsln_RT{2, 2}, 'Correlation \itr');

%% within-subject correlation between pupil baseline and cos angle
% load coeff_robust_pupilBsln_RT
plot_all_data_2tasks(coeff_robust_pupilBsln_cos{1,1}, coeff_robust_pupilBsln_cos{1, 2},...
   coeff_robust_pupilBsln_cos{2, 1}, coeff_robust_pupilBsln_cos{2, 2}, 'Correlation \itr');


%% within-subject correlation between pupil residuals and reaction time 
load coeff_robust_pupilresidphase_RT % group x task
plot_all_data_2tasks(coeff_robust_pupilresidphase_RT{1,1}, coeff_robust_pupilresidphase_RT{1, 2},...
    coeff_robust_pupilresidphase_RT{2, 1}, coeff_robust_pupilresidphase_RT{2, 2}, 'Correlation \itr');

% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(coeff_robust_pupilresidphase_RT{1, 2}, coeff_robust_pupilresidphase_RT{2, 2}, 'Correlation \itr')

% are coefficients different from zero
p_value = cell(2, 2); stats = cell(2, 2);
[h,p_value{1,1},ci,stats{1,1}] = ttest(coeff_robust_pupilresidphase_RT{1,1})
[h,p_value{1,2},ci,stats{1,2}] = ttest(coeff_robust_pupilresidphase_RT{1,2})
[h,p_value{2,1},ci,stats{2,1}] = ttest(coeff_robust_pupilresidphase_RT{2,1})
[h,p_value{2,2},ci,stats{2,2}] = ttest(coeff_robust_pupilresidphase_RT{2,2})

[h,p_value_all,ci,stats_all] = ttest([coeff_robust_pupilresidphase_RT{1,2}, coeff_robust_pupilresidphase_RT{2,2}])

[h,p_value_grp,ci,stats_grp] = ttest2(coeff_robust_pupilresidphase_RT{1,2}, coeff_robust_pupilresidphase_RT{2,2})


%% comparison between correlation coefficients with and without adjustment for ongoing slow fluctuations
load coeff_robust_pupilresidphase_RT % group x task
load coeff_robust_TPR_RT % correlation between TPR amplitude and RT
plot_all_data_2tasks(coeff_robust_TPR_RT{1,1}, coeff_robust_TPR_RT{1, 2},...
    coeff_robust_TPR_RT{2, 1}, coeff_robust_TPR_RT{2, 2}, 'Correlation r');

% difference between coefficients without and with the adjustment
% all participants together
[h,p_value_all,ci,stats_all] = ttest([coeff_robust_pupilresidphase_RT{1,2}, coeff_robust_pupilresidphase_RT{2,2}], ...
    [coeff_robust_TPR_RT{1, 2}, coeff_robust_TPR_RT{2, 2}])



% 
% % are coefficients different across groups
% p_value2 = cell(2, 1); stats2 = cell(2, 1);
% [h,p_value2{1},ci,stats2{1}] = ttest2(coeff_robust_pupilresidphase_RT{1,1}, coeff_robust_pupilresidphase_RT{2,1}); % simple RT
% [h,p_value2{2},ci,stats2{2}] = ttest2(coeff_robust_pupilresidphase_RT{1,2}, coeff_robust_pupilresidphase_RT{2,2}); % gng

% coefficients both groups together one-sample
p_value1 = cell(2, 1); stats1 = cell(2, 1);
[h,p_value1{1},ci,stats1{1}] = ttest([coeff_robust_pupilresidphase_RT{1,1}, coeff_robust_pupilresidphase_RT{2,1}]); % simple RT
[h,p_value1{2},ci,stats1{2}] = ttest([coeff_robust_pupilresidphase_RT{1,2}, coeff_robust_pupilresidphase_RT{2,2}]); % gng

%% within-subject correlation between pupil baseline and TPR amplitude

plot_all_data_2tasks(coeff_robust_pupilBsln_TPR{1,1}, coeff_robust_pupilBsln_TPR{1, 2},...
    coeff_robust_pupilBsln_TPR{2, 1}, coeff_robust_pupilBsln_TPR{2, 2}, 'Correlation r');

%% within-subject correlation between pupil baseline and TPR amplitude residuals
load coeff_robust_pupilBsln_TPRresid
plot_all_data_2tasks(coeff_robust_pupilBsln_TPRresid{1,1}, coeff_robust_pupilBsln_TPRresid{1, 2},...
    coeff_robust_pupilBsln_TPRresid{2, 1}, coeff_robust_pupilBsln_TPRresid{2, 2}, 'Correlation r');

%% within-subject correlation between pupil baseline phase (cos, sin) and RT
cd('G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\PupilVariability');
load coefs_phase_RT % coef1 = cos, coef2 = sin
for coef = 1:3
    plot_all_data_2tasks(coefs_phase_RT{1, 1}(:, coef), coefs_phase_RT{1, 2}(:, coef),...
        coefs_phase_RT{2, 1}(:, coef), coefs_phase_RT{2, 2}(:, coef), 'Coefficient');
    
    % plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
    plot_all_data_onetask(coefs_phase_RT{1, 2}(:, coef), coefs_phase_RT{2, 2}(:, coef), 'Coefficients')
end


[h,p_value,ci,stats] = ttest(coefs_phase_RT{1, 1}(:, 1))
[h,p_value,ci,stats] = ttest(coefs_phase_RT{1, 2}(:, 1))
[h,p_value,ci,stats] = ttest(coefs_phase_RT{2, 1}(:, 1))
[h,p_value,ci,stats] = ttest(coefs_phase_RT{2, 2}(:, 1))

[h,p_value,ci,stats] = ttest(coefs_phase_RT{1, 1}(:, 2))
[h,p_value,ci,stats] = ttest(coefs_phase_RT{1, 2}(:, 2))
[h,p_value,ci,stats] = ttest(coefs_phase_RT{2, 1}(:, 2))
[h,p_value,ci,stats] = ttest(coefs_phase_RT{2, 2}(:, 2))

[h,p_value,ci,stats] = ttest2(coefs_phase_RT{1, 2}(:, 1), coefs_phase_RT{2, 2}(:, 1))
[h,p_value,ci,stats] = ttest2(coefs_phase_RT{1, 2}(:, 2), coefs_phase_RT{2, 2}(:, 2)) 

[h,p_value,ci,stats] = ttest([coefs_phase_RT{1, 2}(:, 1); coefs_phase_RT{2, 2}(:, 1)])
[h,p_value,ci,stats] = ttest([coefs_phase_RT{1, 2}(:, 2); coefs_phase_RT{2, 2}(:, 2)])


%% within-subject correlation between pupil baseline phase (cos, sin) and RT controlling for time-on-task
load coefs_phase_timeontask_RT
for coef = 1:4
    plot_all_data_2tasks(coefs_phase_timeontask_RT{1, 1}(:, coef), coefs_phase_timeontask_RT{1, 2}(:, coef),...
        coefs_phase_timeontask_RT{2, 1}(:, coef), coefs_phase_timeontask_RT{2, 2}(:, coef), 'Coefficient');
end

[h,p_value,ci,stats] = ttest([coefs_phase_timeontask_RT{1, 1}(:, 1); coefs_phase_timeontask_RT{2, 1}(:, 1)])

[h,p_value,ci,stats] = ttest([coefs_phase_timeontask_RT{1, 2}(:, 1); coefs_phase_timeontask_RT{2, 2}(:, 1)])


%% variability of PD responses after regressing out the effect of ongoing signal fluctuations
load pupil_resd_phase_std

plot_all_data_2tasks(pupil_resd_phase_std{1,1}', pupil_resd_phase_std{1, 2}',...
    pupil_resd_phase_std{2, 1}', pupil_resd_phase_std{2, 2}', 'PD residuals SD');

% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(pupil_resd_phase_std{1, 2}', pupil_resd_phase_std{2, 2}', 'PD residuals SD')

[h,p_value,ci,stats] = ttest2(pupil_resd_phase_std{1, 2}', pupil_resd_phase_std{2, 2}')


%% are coefficients different from zero
p_value = cell(2, 2); stats = cell(2, 2);
[h,p_value{1,1},ci,stats{1,1}] = ttest(coeff_robust_pupilBsln_TPR{1,1})
[h,p_value{1,2},ci,stats{1,2}] = ttest(coeff_robust_pupilBsln_TPR{1,2})
[h,p_value{2,1},ci,stats{2,1}] = ttest(ccoeff_robust_pupilBsln_TPR{2,1})
[h,p_value{2,2},ci,stats{2,2}] = ttest(coeff_robust_pupilBsln_TPR{2,2})

% are coefficients different across groups
p_value2 = cell(2, 1); stats2 = cell(2, 1);
[h,p_value2{1},ci,stats2{1}] = ttest2(coeff_robust_pupilBsln_TPR{1,1}, coeff_robust_pupilBsln_TPR{2,1}); % simple RT
[h,p_value2{2},ci,stats2{2}] = ttest2(coeff_robust_pupilBsln_TPR{1,2}, coeff_robust_pupilBsln_TPR{2,2}); % gng

% coefficients both groups together one-sample
p_value1 = cell(2, 1); stats1 = cell(2, 1);
[h,p_value1{1},ci,stats1{1}] = ttest([coeff_robust_pupilBsln_TPR{1,1}, coeff_robust_pupilBsln_TPR{2,1}]); % simple RT
[h,p_value1{2},ci,stats1{2}] = ttest([coeff_robust_pupilBsln_TPR{1,2}, coeff_robust_pupilBsln_TPR{2,2}]); % gng






%% create excel file for SPSS analysis
load coeff_robust_TPR_RT; load coeff_robust_pupilBsln_TPR;
load coefs_phase_RT; load coeff_robust_pupilresidphase_RT;
load coeff_robust_pupilBsln_TPR; load pupil_resd_phase_std
clear T_values
T_values(:, 1) = [ones(length(coeff_robust_pupilresidphase_RT{1, 1}), 1); ones(length(coeff_robust_pupilresidphase_RT{2, 1}), 1)*2];
T_values =  [T_values, [coeff_robust_pupilresidphase_RT{1, 1}'; coeff_robust_pupilresidphase_RT{2, 1}']];
T_values =  [T_values, [coeff_robust_pupilresidphase_RT{1, 2}'; coeff_robust_pupilresidphase_RT{2, 2}']];
T_values =  [T_values, [coeff_robust_pupilBsln_TPR{1, 1}'; coeff_robust_pupilBsln_TPR{2, 1}']];
T_values =  [T_values, [coeff_robust_pupilBsln_TPR{1, 2}'; coeff_robust_pupilBsln_TPR{2, 2}']];

T_values =  [T_values, [pupil_resd_phase_std{1, 1}'; pupil_resd_phase_std{2, 1}']];
T_values =  [T_values, [pupil_resd_phase_std{1, 2}'; pupil_resd_phase_std{2, 2}']];

% estimated coefs between baseline phase and RT
T_values =  [T_values, [coefs_phase_RT{1, 1}(:, 1:2); coefs_phase_RT{2, 1}(:, 1:2)]];
T_values =  [T_values, [coefs_phase_RT{1, 2}(:, 1:2); coefs_phase_RT{2, 2}(:, 1:2)]];

% correlation coefficients between PD amplitude and RT
T_values =  [T_values, [coeff_robust_TPR_RT{1, 1}'; coeff_robust_TPR_RT{2, 1}']];
T_values =  [T_values, [coeff_robust_TPR_RT{1, 2}'; coeff_robust_TPR_RT{2, 2}']];

all_column_names = {'group','corr_pupilres_simpleRT','corr_pupilres_gng',...
    'corr_pupilbsln_simpleRT','corr_pupilbsln_gng', ...
    'pupilres_std_simpleRT', 'pupilres_std_gng', ...
    'cos_RT_simpleRT', 'sin_RT_simpleRT', 'cos_RT_gng', 'sin_RT_gng', ...
    'corr_PD_RT_SimpleRT', 'corr_PD_RT_gng'};

T = array2table(T_values, ...
    'VariableNames',all_column_names);

filename = 'RTcorr_PupilBsln_PupilResid.xlsx';
writetable(T,filename,'Sheet',1,'Range','A1')



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
    
    % plot line on zero
    plot([0 6], zeros(2, 1), '--k')
    

    % axes('XColor','none');
    hold off;
    axis([0 6 -inf inf]);
    ax = gca;
    c = ax.Color;
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel=[1 2 1 2];
    xticks([1 2 4 5])
    ylabel(y_label_text, 'FontSize',32, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end

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
   

    % group 2
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
    ax.FontSize = 18;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel= [];
%     xticks([1 2])
    ylabel(y_label_text, 'FontSize', 24, 'FontWeight','normal')
    title(title_text, 'FontSize', 24, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end


function plot_data_mean_std_young_older(data_young, data_older, graph_title, y_axis_label)

mean_data_young=mean(data_young, 1);
se_data_young=std(data_young, 0, 1)/sqrt(size(data_young, 1));

mean_data_older=mean(data_older, 1);
se_data_older=std(data_older, 0, 1)/sqrt(size(data_older, 1));


% color gng young = [0 .5 0]
% color simple RT young =   [1 .5 0]
% color gng older = [0 0 .5]
% color simple RT older =  [.75 0.25 0]

x_axis=-.999:1/240:6;
figure;
for x = 1.001:1/240:1.500 % grey background between 1000 amd 1500ms after ceu-onset
    plot([x x],[-.009 .099], 'color', [.8 .8 .8] ); hold on
end
hold on
% detection data
plot(x_axis, mean_data_young,'color', 'k', 'LineWidth',1.5)
hold on
jbfill(x_axis,mean_data_young+se_data_young, mean_data_young-se_data_young, 'k','k', 0.1)
% GNG data
hold on
plot(x_axis, mean_data_older,'color', 'r', 'LineWidth',1.5)
hold on
jbfill(x_axis, mean_data_older+se_data_older, mean_data_older-se_data_older,'r','r', 0.1)
hold on
plot(x_axis, zeros(length(x_axis), 1), '--k')
hold off
ax = gca;
% c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 24;
ax.FontName = 'Arial';
ax.Color = 'none';
axis([-.2 4 -.01 .1])
xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
ylabel(y_axis_label, 'FontSize', 32, 'FontWeight','normal')
title(graph_title, 'FontSize', 32, 'FontWeight','normal')

end


function plot_quintiles(data, title_text)

    Mean1= squeeze(mean(data(:,1,:), 1));
    Mean2= squeeze(mean(data(:,2,:), 1));
    Mean3= squeeze(mean(data(:,3,:), 1));
    Mean4= squeeze(mean(data(:,4,:), 1));
    Mean5= squeeze(mean(data(:,5,:), 1));

    SE1=(squeeze(std(data(:,1,:), 0, 1))/sqrt(size(data, 1)));
    SE2=(squeeze(std(data(:,2,:), 0, 1))/sqrt(size(data, 1)));
    SE3=(squeeze(std(data(:,3,:), 0, 1))/sqrt(size(data, 1)));
    SE4=(squeeze(std(data(:,4,:), 0, 1))/sqrt(size(data, 1)));
    SE5=(squeeze(std(data(:,5,:), 0, 1))/sqrt(size(data, 1)));

    colormap cool;%colormap summer;
    cmap = colormap;

    figure;
    xaxis = -.999:1/240:6;
    % zero line
    plot(xaxis, zeros(length(xaxis)), 'k:');
    hold on
    plot( xaxis', Mean1, 'color', cmap(10,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean1+SE1)',(Mean1-SE1)', cmap(10,:), cmap(10,:), 1, 0.2)
    hold on
    plot( xaxis, Mean2, 'color', cmap(20,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean2+SE2)',(Mean2-SE2)', cmap(20,:), cmap(20,:), 1, 0.2)
    hold on
    plot( xaxis, Mean3, 'color', cmap(30,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean3+SE3)',(Mean3-SE3)', cmap(30,:), cmap(30,:), 1, 0.2)
    hold on
    plot( xaxis, Mean4, 'color', cmap(40,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean4+SE4)',(Mean4-SE4)', cmap(40,:), cmap(40,:), 1, 0.2)
    hold on
    plot( xaxis, Mean5, 'color', cmap(50,:), 'linewidth', 2);
    hold on
    jbfill(xaxis,(Mean5+SE5)',(Mean5-SE5)', cmap(50,:), cmap(50,:), 1, 0.2)
    hold on
    % plot(xaxis, zeros(1, 840), ':k')
    hold off
    axis([0 2 -.005 .08]);
    % title('GNG', 'FontSize', 24, 'FontWeight','bold')
    % legend('RT1', 'RT2', 'RT3', 'RT4', 'RT5', 'location', 'northwest')
    ax = gca;
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Pupil response (%)', 'FontSize', 32, 'FontWeight','normal')
    title(title_text, 'FontSize', 32, 'FontWeight','normal');
    
end



function plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)

    % plot data for young group - go/nogo task
    figure; box off; hold on
    plot([0 3], [0 0], '--k'); % line at zero
    hold on
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
    axis([0 3 -inf inf]);
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

