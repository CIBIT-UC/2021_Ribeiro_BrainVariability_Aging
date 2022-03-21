% adjusting ERP amplitude for phase of ongoing slow fluctuations and
% checking correlation with RT
clear; close all;
% participants id
% eeg included participants
young_eeg=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older_eeg=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];
% pupil included participants
young_pupil = [4	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46 50 51 53 54 56	59	62 66	68	72	74	76  78  80 81   82  84  85];
older_pupil = [7  8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77 79   83  86];
% find participants with eeg and pupil data
young = mintersect(young_eeg, young_pupil);
older = mintersect(older_eeg, older_pupil);
group = {young, older};
task={'D1', 'D2', 'G1', 'G2'};

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

coeff_robust_chan = cell(2, 2);
for grp = 1:2
    part = 0;
    for p = group{grp}
        part = part + 1;
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
             
        % erp data
        eeg_directory = strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\EEG\');
        epochs_erp = {}; RT_eeg = {}; RT_pupil = {}; ERPs = {}; EEG_phase_bsln = {}; EEG_amplit_bsln = {};
        
       for t = 1:4 % task runs 2 simple RT + 2 gng
             % clear eeglab
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            %load file - from eeglab_analysis_2b_synchronize_pupil
            % _RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica_ArtICsRem_Filt0_1_35Hz_Pupil.set
            % manual/visual removal of periods with artifacts
            % save as ..._Filt0_1_35Hz_Pupil_ManualArtRej.set
            if p == 77 && ismember(t, [1, 2])
                filename=strcat('AB', num2str(p), '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej_HEO.set');
            else
                filename=strcat('AB', num2str(p), '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej.set');
            end
            EEG = pop_loadset(filename, eeg_directory);
            % epoch data
            EEG = pop_epoch( EEG, {  '1'  }, [-.2  6]);
           
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
            % with urevent - note that error trials are very rare
            % error_trials - all errors, response2cue, multiple responses, misses, slow responses
            error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses, response2nogo);
            % only correct go trials, excluding trials after error
            include_epochs = setdiff(correct_trial, error_trials+1);
            
            cnt = 0;
            for inc = include_epochs'
                cnt = cnt + 1;
                RT_eeg{t}(cnt, 1) = RT(RT(:, 1) == inc, 2);
    %             cue2target_interval_inc{t}(cnt, 1) = cue2target_interval(cue2target_interval(:, 1) == inc, 2);
            end
            
            
            EEG = pop_select( EEG, 'trial',include_epochs);
%             % create variable with average EEG amplitude 200 ms baseline
%             % before cue onset
%             EEGbsln{t} = squeeze(mean(EEG.data(1:59, 1:99, :), 2)); 
            
            EEG = pop_rmbase( EEG, [-200 0] ,[]);
           % ERP data - avg amplitude between 1000 and 1500 ms after cue onset
           ERPs{t} = squeeze(mean(EEG.data(1:59, 600:850, :), 2)); % chan X trials    
%            ERPs_single_trial{t} = EEG.data(1:59, 1:1100, :); % chan x time x trial
           
           % calculate phase of slow fluctuations at cue onset
           EEG = pop_loadset(filename, eeg_directory);
            % interpolate bad channels removed earlier
            if length(EEG.chanlocs)<67 && p ~= 6 && p~=38 % these participants do not have eyetracking data but have all channels ok
                [EEG] = interpol( EEG);
            end
            % filter data low pass 4 Hz
            EEG = pop_eegfiltnew(EEG, 'hicutoff',4,'plotfreqz',0);
            % hilbert transform
            signal_phase = [];
            for chan = 1:59
                y = hilbert(EEG.data(chan, :));
                signal_phase(chan, :) = angle(y); % extract instantaneous phase
            end

            EEG.data(1:59, :) = signal_phase;
            % epoch data
            EEG = pop_epoch( EEG, {  '1'  }, [-.2  6]);
            EEG = pop_select( EEG, 'trial',include_epochs); % exclude errors etc
           % phase at cue-onset
           EEG_phase_bsln{t} = squeeze(EEG.data(1:59, 100, :)); % chan x trials
           
           
           % calculate amplitude envelop of slow fluctuations at cue onset
           EEG = pop_loadset(filename, eeg_directory);
            % interpolate bad channels removed earlier
            if length(EEG.chanlocs)<67 && p ~= 6 && p~=38 % these participants do not have eyetracking data but have all channels ok
                [EEG] = interpol( EEG);
            end
            % filter data low pass 4 Hz
            EEG = pop_eegfiltnew(EEG, 'hicutoff',4,'plotfreqz',0);
            % hilbert transform
            signal_amplitude = [];
            for chan = 1:59
                y = hilbert(EEG.data(chan, :));
                signal_amplitude(chan, :) = abs(y); % extract instantaneous amplitude envelop
            end

            EEG.data(1:59, :) = signal_amplitude;
            % epoch data
            EEG = pop_epoch( EEG, {  '1'  }, [-.2  6]);
            EEG = pop_select( EEG, 'trial', include_epochs); % exclude errors etc
           % amplitude envelope at cue-onset
           EEG_amplit_bsln{t} = squeeze(EEG.data(1:59, 100, :)); % chan x trials
           
           % determine original epoch number to match with pupil data
           % for each original epoch find the event original event (cue) number 
            a=0; trial_numbers=[]; trial_numbers(:, 1)=1:60; % there are 60 trials in each run - original trial number
            for i=1:size(EEG.urevent, 2)
                if EEG.urevent(i).type == 1
                    a=a+1;
                    trial_numbers(a, 2) = i; % event number on original file without periods with artefacts removed
                elseif strcmp(EEG.urevent(i).type, {'boundary'})
                    trial_numbers_boundary = i;
                end
            end
           epochs_erp{t} = zeros(size(EEG.epoch, 2), 2); % epoch number x epoch original
           epochs_erp{t}(:, 1) = 1:size(EEG.epoch, 2);
           for epoch = 1:size(EEG.epoch, 2)
                for i = 1:size(EEG.event, 2)
                    if EEG.event(i).epoch == epoch && strcmp(EEG.event(i).type, '1')
                        epochs_erp{t}(epoch, 2) = trial_numbers((trial_numbers(:, 2) == EEG.event(i).urevent), 1);
                    end
                end
           end
        
       end 
       
       
       %% calculate pupil evoked response
        pupil_directory=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\');
        load([pupil_directory, 'RejectEpochsCorrect']); % data from G:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\pupil_analysis_3b_RejectedEpochs_SaveRunData.m
        
        clear pupil_avg_response_tmp pupil_response
        RT_pupil = cell(4, 1); time_on_task = cell(4, 1);
        pupil_avg_response = cell(4, 1); pupil_bsln_phase = cell(4, 1); pupil_bsln_amp = cell(4, 1);
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
                    RT_pupil{t}=[RT_pupil{t}; (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/240]; 
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
            
       %% match epochs of pupil, alpha, and erp
%        pupil_epochs
%        epochs_prestim_spectra
%        epochs_erp
       included_epochs = {}; pupil_epochs = time_on_task;
       for t = 1:4
            included_epochs{t} = mintersect(pupil_epochs{t}, epochs_erp{t}(:, 2));

            % delete epochs not to include
            % pupil and reaction time
            index2delete = []; epochs2delete = setdiff(pupil_epochs{t}, included_epochs{t});
            if ~isempty(epochs2delete)
                for epoch = 1:length(epochs2delete)
                    index2delete = [index2delete; find(pupil_epochs{t} == epochs2delete(epoch))];
                end
                pupil_avg_response{t}(index2delete) = [];
                pupil_bsln_phase{t}(index2delete) = [];
                pupil_bsln_amp{t}(index2delete) = [];
                RT_pupil{t}(index2delete) = [];
            end

            % erp
            index2delete = []; epochs2delete = setdiff(epochs_erp{t}(:, 2), included_epochs{t});
            if ~isempty(epochs2delete)
                for epoch = 1:length(epochs2delete)
                    index2delete = [index2delete; find(epochs_erp{t}(:, 2) == epochs2delete(epoch))];
                end
                ERPs{t}(:, index2delete) = [];
                EEG_amplit_bsln{t}(:, index2delete) = [];
                EEG_phase_bsln{t}(:, index2delete) = [];
                RT_eeg{t}(index2delete) = [];
            end
       end
       
      %% adjust ERP amplitude by reggressing out phase of slow fluctuations at baseline
      % tasks: simple RT (runs 1 and 2) and gng (runs 3 and 4) 
       Res_erp_simpleRT = []; Res_erp_gng = [];
       for chan  = 1:59
           % regress eeg phase effect at baseline out of erp amplitude
           % simple RT
           y = [ERPs{1}(chan, :)'; ERPs{2}(chan, :)'];
           x = [[EEG_amplit_bsln{1}(chan, :)'; EEG_amplit_bsln{2}(chan, :)'].*cos([EEG_phase_bsln{1}(chan, :)'; EEG_phase_bsln{2}(chan, :)']),...
               [EEG_amplit_bsln{1}(chan, :)'; EEG_amplit_bsln{2}(chan, :)'].*sin([EEG_phase_bsln{1}(chan, :)'; EEG_phase_bsln{2}(chan, :)']), ones(length(y), 1)];
           [~,~,Res_erp_simpleRT(:, chan)] = regress(y, x);
           
           % gng
           y = [ERPs{3}(chan, :)'; ERPs{4}(chan, :)'];
           x = [[EEG_amplit_bsln{3}(chan, :)'; EEG_amplit_bsln{4}(chan, :)'].*cos([EEG_phase_bsln{3}(chan, :)'; EEG_phase_bsln{4}(chan, :)']),...
               [EEG_amplit_bsln{3}(chan, :)'; EEG_amplit_bsln{4}(chan, :)'].*sin([EEG_phase_bsln{3}(chan, :)'; EEG_phase_bsln{4}(chan, :)']), ones(length(y), 1)];
           [~,~,Res_erp_gng(:, chan)] = regress(y, x);
       end
       
       % regress phase of slow fluctuations out of pupil response amplitude
       % simple RT
       y = [pupil_avg_response{1}; pupil_avg_response{2}];
       x = [[pupil_bsln_amp{1}.*cos(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*cos(pupil_bsln_phase{2})], [pupil_bsln_amp{1}.*sin(pupil_bsln_phase{1}); pupil_bsln_amp{2}.*sin(pupil_bsln_phase{2})], ones(length(y), 1)];
       [~,~,Res_pupil_simpleRT] = regress(y, x);

       % gng
       y = [pupil_avg_response{3}; pupil_avg_response{4}];
       x = [[pupil_bsln_amp{3}.*cos(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*cos(pupil_bsln_phase{4})], [pupil_bsln_amp{3}.*sin(pupil_bsln_phase{3}); pupil_bsln_amp{4}.*sin(pupil_bsln_phase{4})], ones(length(y), 1)];
       [~,~,Res_pupil_gng] = regress(y, x);
       
       % compare RT from eeg with RT from pupil data
       RT_diff{grp}{part} = [[RT_eeg{1}; RT_eeg{2}; RT_eeg{3}; RT_eeg{4}], [RT_pupil{1}; RT_pupil{2}; RT_pupil{3}; RT_pupil{4}], [RT_eeg{1}; RT_eeg{2}; RT_eeg{3}; RT_eeg{4}]-[RT_pupil{1}; RT_pupil{2}; RT_pupil{3}; RT_pupil{4}]];
       
%        % correlate erp and pupil response before adjustment
%        for chan = 1:size(ERPs{1}, 1)
%            % simple RT
%             [r,~,~,~,~,~] = skipped_correlation([ERPs{1}(chan, :)'; ERPs{2}(chan, :)'], [pupil_avg_response{1}; pupil_avg_response{2}],0);
%             coeff_erp_pupil{grp, 1}(part, chan) = r.Pearson;
%              % gng
%             [r,~,~,~,~,~] = skipped_correlation([ERPs{3}(chan, :)'; ERPs{4}(chan, :)'], [pupil_avg_response{3}; pupil_avg_response{4}],0);
%             coeff_erp_pupil{grp, 2}(part, chan) = r.Pearson;
%        end
%        
       
       %%
         % create variables for multiple regression
         % ERP residuals for each channel and pupil response residuals as
         % independent variables
         % dependent variable RT
         % simple RT
        y = [RT_eeg{1}; RT_eeg{2}];
        for chan = 1:size(ERPs{1}, 1)
            x = [Res_erp_simpleRT(:, chan), Res_pupil_simpleRT, ones(length(y), 1)];
%           [b,bint,r,rint,stats] = regress(y,X)
            [coeff_rt_erp_pupil_phaseresid{grp, 1}(part, chan, :),~,~,~,stats_rt_erp_pupil_phaseresid{grp, 1}(part, chan, :)] = regress(y,x);
            [coeff_rt_erp_phaseresid{grp, 1}(part, chan, :),~,~,~,stats_rt_erp_phaseresid{grp, 1}(part, chan, :)] = regress(y,x(:, [1 3]));
            % within-subject correlation between pupil response amplitude
            % residuals and ERP residuals
%             [r,~,~,~,~,~] = skipped_correlation(Res_erp_simpleRT(:, chan), Res_pupil_simpleRT,0);
%             coeff_erp_pupil_phaseresid{grp, 1}(part, chan) = r.Pearson;
        end
        [coeff_rt_pupil_phaseresid{grp, 1}(part, :),~,~,~,stats_rt_pupil_phaseresid{grp, 1}(part, :)] = regress(y,x(:, [2 3]));

        % gng
        y = [RT_eeg{3}; RT_eeg{4}];
        for chan = 1:size(ERPs{1}, 1)
            x = [ Res_erp_gng(:, chan), Res_pupil_gng, ones(length(y), 1)];
            [coeff_rt_erp_pupil_phaseresid{grp, 2}(part, chan, :),~,~,~,stats_rt_erp_pupil_phaseresid{grp, 2}(part, chan, :)] = regress(y,x);
            [coeff_rt_erp_phaseresid{grp, 2}(part, chan, :),~,~,~,stats_rt_erp_phaseresid{grp, 2}(part, chan, :)] = regress(y,x(:, [1 3]));
            % within-subject correlation between pupil response amplitude
            % residuals and ERP residuals
%             [r,~,~,~,~,~] = skipped_correlation(Res_erp_gng(:, chan), Res_pupil_gng,0);
%             coeff_erp_pupil_phaseresid{grp, 2}(part, chan) = r.Pearson;
            if chan == 25
                figure; plot(Res_erp_gng(:, chan), Res_pupil_gng, 'o');
            end
        end  
        [coeff_rt_pupil_phaseresid{grp, 2}(part, :),~,~,~,stats_rt_pupil_phaseresid{grp, 2}(part, :)] = regress(y,x(:, [2 3]));
    end
end
save coeff_rt_erp_pupil_phaseresid coeff_rt_erp_pupil_phaseresid
save stats_rt_erp_pupil_phaseresid stats_rt_erp_pupil_phaseresid
save coeff_rt_erp_phaseresid coeff_rt_erp_phaseresid
save stats_rt_erp_phaseresid stats_rt_erp_phaseresid
save coeff_rt_pupil_phaseresid coeff_rt_pupil_phaseresid
save stats_rt_pupil_phaseresid stats_rt_pupil_phaseresid

% save coeff_erp_pupil_phaseresid coeff_erp_pupil_phaseresid
% save coeff_erp_pupil coeff_erp_pupil

%% check RT differences = sanity check = all ok
% RT_diff{grp, part}
figure; 
for grp = 1:2
    for part = 1:length(RT_diff{grp})
        plot(RT_diff{grp}{part}(:, 1), RT_diff{grp}{part}(:, 2), 'o')
        waitforbuttonpress
    end
end

%% multiple regression rt = erp and pupil response signals
load coeff_rt_erp_pupil_phaseresid % grp task - part, chan, coeffs
task_name = {'Simple RT' 'Go/no-go'};

clear pval_t1 t_orig_t1 crit_t_t1 pval_t1_all t_orig_t1_all crit_t_t1_all
for task = 1:2
    for coef = 1:2
        for grp = 1:2
            [pval_t1(grp, task, coef, :), t_orig_t1(grp, task, coef, :), crit_t_t1(grp, task, coef, :),~,~]=mult_comp_perm_t1(squeeze(coeff_rt_erp_pupil_phaseresid{grp, task}(:, :, coef)));
        end
        [pval_t2(task, :, coef), t_orig_t2(task, :, coef), crit_t_t2(task, :, coef),~,~] = mult_comp_perm_t2(squeeze(coeff_rt_erp_pupil_phaseresid{1, task}(:, :, coef)),squeeze(coeff_rt_erp_pupil_phaseresid{2, task}(:, :, coef)));
        [pval_t1_all(task, coef, :), t_orig_t1_all(task, coef, :), crit_t_t1_all(task, coef, :),~,~]=mult_comp_perm_t1([squeeze(coeff_rt_erp_pupil_phaseresid{1, task}(:, :, coef)); squeeze(coeff_rt_erp_pupil_phaseresid{2, task}(:, :, coef))]);
        
        plot_data_all_electrodes(squeeze(coeff_rt_erp_pupil_phaseresid{1, task}(:, :, coef)), squeeze(coeff_rt_erp_pupil_phaseresid{2, task}(:, :, coef)), 'Regression Coeffs', ['Coef', num2str(coef), ' - ', task_name(task)])
    end
end

% plot t-values from one sample t-test both groups together
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_txt = {'Simple RT' 'Go/no-go'};%cmap = crameri('batlow');
for task = 1:2
    for coef = 1:2
        sig_t_values = zeros(1, 59);
        sig_t_values(t_orig_t1_all(task, coef, :) > crit_t_t1_all(task, coef, 1)) = 1;
        sig_t_values(t_orig_t1_all(task, coef, :) < crit_t_t1_all(task, coef, 2)) = 1;
        sig_chan_number = find(sig_t_values == 1);

        for x = 1:length(sig_chan_number)
            sig_chans_all{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
        end

        figure; 
        if ~isempty(sig_chan_number)
            topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
        end
        topoplot(t_orig_t1_all(task, coef, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        if coef == 1
            caxis([0 5]); %
        else
            caxis([-9 0]);
        end
        colorbar; c.Axis.FontSize = 16;
%         colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
%         colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
        title(title_txt{task}, 'FontSize', 30, 'FontWeight','normal')
        set(get(gca,'title'),'Position',[0,-.65, 0])


% %         plot average correlation coefficients with significant electrodes
% %         highlighted
% %         figure; 
% %         if ~isempty(sig_chan_number)
% %             topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
% %                     'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
% %         end
% %         topoplot(mean([coeff_rt_erp_pupil_phaseresid{1, task}(:, :, coef); coeff_rt_erp_pupil_phaseresid{2, task}(:, :, coef)], 1), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
% %         caxis([0 .1]); %c.Axis.FontSize = 16;
% %         colorbar;
% %         colorbar('Ticks',[0, .05, .1], 'FontSize', 30, 'FontWeight','normal');
% %         colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
    end
end


%% check how much R2 of models change with including pupil in them
load stats_rt_pupil_phaseresid; load stats_rt_erp_phaseresid; load stats_rt_erp_pupil_phaseresid
for grp = 1:2
    for task = 1:2
        R2_pupil{grp, task} = squeeze(stats_rt_pupil_phaseresid{grp, task}(:, 1));
        R2_erp{grp, task} = squeeze(stats_rt_erp_phaseresid{grp, task}(:, :, 1));
        R2_erp_pupil{grp, task} = squeeze(stats_rt_erp_pupil_phaseresid{grp, task}(:, :, 1));
    end
end

for task = 1:2
    for grp = 1:2
        [pval_effect_of_pupil(grp, task, :), t_orig_effect_of_pupil(grp, task, :), crit_t_effect_of_pupil(grp, task, :),~,~]=mult_comp_perm_t2(R2_erp_pupil{grp, task}, R2_erp{grp, task});
        [pval_effect_of_erp(grp, task, :), t_orig_effect_of_erp(grp, task, :), crit_t_effect_of_erp(grp, task, :),~,~]=mult_comp_perm_t2(R2_erp_pupil{grp, task}, repmat(R2_pupil{grp, task}, 1, 59));
    end
end

% plot t-values from t-test
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_txt = {'Simple RT' 'Go/no-go'}; grp_txt = {'young', 'older'};
for ef = 1:2 %pupil , erp
    for grp =1:2
        for task = 2%1:2
            sig_t_values = zeros(1, 59);
            if ef == 1
                sig_t_values(t_orig_effect_of_pupil(grp, task, :) > crit_t_effect_of_pupil(grp, task, 2)) = 1;
                sig_t_values(t_orig_effect_of_pupil(grp, task, :) < crit_t_effect_of_pupil(grp, task, 1)) = 1;
            else
                sig_t_values(t_orig_effect_of_erp(grp, task, :) > crit_t_effect_of_erp(grp, task, 2)) = 1;
                sig_t_values(t_orig_effect_of_erp(grp, task, :) < crit_t_effect_of_erp(grp, task, 1)) = 1; 
            end
            sig_chan_number = find(sig_t_values == 1);

            for x = 1:length(sig_chan_number)
                sig_chans_all{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
            end

            figure; 
            if ~isempty(sig_chan_number)
                topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                        'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
            end
            if ef == 1
                topoplot(t_orig_effect_of_pupil(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
            else
                topoplot(t_orig_effect_of_erp(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off');      
            end
        %     if coef == 1
        %         caxis([0 5]); %
        %     else
        %         caxis([-9 0]);
        %     end
            colorbar; c.Axis.FontSize = 16;
%                 colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
            colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
            title([title_txt{task}, grp_txt{grp}], 'FontSize', 30, 'FontWeight','normal')
            set(get(gca,'title'),'Position',[0,-.65, 0])
        end
    end
end

%% compare R2 at channel C5 with and without pupil
% matrix with R2 pupil, R2 erp, R2 pupil + erp - gng only
% find chan index = C5
for chan = 1:length(chanlocs_EEGChanOnly)
    if strcmp(chanlocs_EEGChanOnly(chan).labels, 'C5')
        c5 = chan;
    end
end

R2_C5 = [[ones(length(R2_pupil{1, 2}), 1); ones(length(R2_pupil{2, 2}), 1)*2], ...
    [R2_pupil{1, 2}; R2_pupil{2, 2}], [R2_erp{1, 2}(:, c5); R2_erp{2, 2}(:, c5)],...
    [R2_erp_pupil{1, 2}(:, c5); R2_erp_pupil{2, 2}(:, c5)]];

% find channel with strongest R2 in model only with ERP
[B_all, I_sort_all] = sortrows(mean([mean(R2_erp{1, 2}, 1)', mean(R2_erp{2, 2}, 1)'], 2));
chanlocs_EEGChanOnly(I_sort_all(end)).labels
% I_sort_all(end) = FC6
R2_FC6 = [[R2_erp{1, 2}(:, I_sort_all(end)); R2_erp{2, 2}(:, I_sort_all(end))],...
    [R2_erp_pupil{1, 2}(:, I_sort_all(end)); R2_erp_pupil{2, 2}(:, I_sort_all(end))]];

R2_avg = [[mean(R2_erp{1, 2}, 2); mean(R2_erp{2, 2}, 2)],...
    [mean(R2_erp_pupil{1, 2}, 2); mean(R2_erp_pupil{2, 2}, 2)]];

R2 = array2table([R2_C5, R2_FC6, R2_avg], "VariableNames",["group","pupil","erp_c5", "pupil_erp_c5", "erp_fc6", "pupil_erp_fc6",  "erp_avg", "pupil_erp_avg"]);

writetable(R2, 'R2_adjusted_pupil_erp_vs_RT.xlsx')




%% check if association between RT and ERP is significant
load coeff_rt_erp_phaseresid % coeff_rt_erp_phaseresid{grp, 2}(part, chan, :)
task_name = {'Simple RT' 'Go/no-go'};

clear pval_t1 t_orig_t1 crit_t_t1 pval_t1_all t_orig_t1_all crit_t_t1_all pval_t2 t_orig_t2 crit_t_t2
for task = 1:2
    for grp = 1:2
        [pval_t1(grp, task, :), t_orig_t1(grp, task, :), crit_t_t1(grp, task, :),~,~]=mult_comp_perm_t1(squeeze(coeff_rt_erp_phaseresid{grp, task}(:, :, 1)));
    end
    [pval_t2(task, :), t_orig_t2(task, :), crit_t_t2(task, :),~,~] = mult_comp_perm_t2(squeeze(coeff_rt_erp_phaseresid{1, task}(:, :, 1)),squeeze(coeff_rt_erp_phaseresid{2, task}(:, :, 1)));
    [pval_t1_all(task, :), t_orig_t1_all(task, :), crit_t_t1_all(task, :),~,~]=mult_comp_perm_t1([squeeze(coeff_rt_erp_phaseresid{1, task}(:, :, 1)); squeeze(coeff_rt_erp_phaseresid{2, task}(:, :, 1))]);
end

% plot t-values from one sample t-test both groups together
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_txt = {'Simple RT' 'Go/no-go'};%cmap = crameri('batlow');
for task = 1:2
    sig_t_values = zeros(1, 59);
    sig_t_values(t_orig_t1_all(task, :) > crit_t_t1_all(task, 1)) = 1;
    sig_t_values(t_orig_t1_all(task, :) < crit_t_t1_all(task, 2)) = 1;
    sig_chan_number = find(sig_t_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans_all{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end

    figure; 
    if ~isempty(sig_chan_number)
        topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
    end
    topoplot(t_orig_t1_all(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
%     if coef == 1
%         caxis([0 5]); %
%     else
%         caxis([-9 0]);
%     end
    colorbar; c.Axis.FontSize = 16;
%     colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
    title(title_txt{task}, 'FontSize', 30, 'FontWeight','normal')
    set(get(gca,'title'),'Position',[0,-.65, 0])
end



%% correlation between erp and pupil response with and without adjustment
title_txt = {'Simple RT' 'Go/no-go'}; adj_text = {'with adj' 'without adj'}; grp_txt = {'young', 'older'};
clear pval_t2 t_orig_t2 crit_t_t2 sig_chans_all
for adj = 1:2
    if adj == 1 % with adjustment
        load coeff_erp_pupil_phaseresid; % group x task
        coeffs = coeff_erp_pupil_phaseresid;
    else % without adjustment
        load coeff_erp_pupil;
        coeffs = coeff_erp_pupil;
    end
    clear pval_t1_all t_orig_t1_all crit_t_t1_all pval_t1 t_orig_t1 crit_t_t1
    clear pval_t1 t_orig_t1 crit_t_t1
    for task = 1:2
        for grp = 1:2
            [pval_t1(grp, task, :), t_orig_t1(grp, task, :), crit_t_t1(grp, task, :),~,~]=mult_comp_perm_t1(coeffs{grp, task});
        end
    end
    
    % plot t-values from one sample t-test 
    load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
    title_txt = {'Simple RT' 'Go/no-go'};%cmap = crameri('batlow');
    for grp = 1:2
        for task = []%2%1:2

            sig_t_values = zeros(1, 59);
            sig_t_values(t_orig_t1(grp, task, :) > crit_t_t1(grp, task, 2)) = 1;
            sig_t_values(t_orig_t1(grp, task, :) < crit_t_t1(grp, task, 1)) = 1;
            sig_chan_number = find(sig_t_values == 1);

            for x = 1:length(sig_chan_number)
                sig_chans_all{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
            end

            figure;
            if ~isempty(sig_chan_number)
                topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                        'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
            end
            topoplot(mean(coeffs{grp, task}, 1), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
            caxis([-.18 0]); %c.Axis.FontSize = 16;
            colorbar;
            colorbar('Ticks',[-.15, -.1, -.05, 0], 'FontSize', 26, 'FontWeight','normal');
            colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
%             title([grp_txt{grp}, '-',adj_text{adj}], 'FontSize', 30, 'FontWeight','normal')
%             set(get(gca,'title'),'Position',[0,-.65, 0])
        end
    end
    
    % are the correlation coefficients different across groups?
    for task = 1:2
        [pval_t2(adj, task, :), t_orig_t2(adj, task, :), crit_t_t2(adj, task, :),~,~] = mult_comp_perm_t2(coeffs{1, task}, coeffs{2, task});
    end
    
     % plot t-values
    load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat

    for task = 2%1:2
        sig_t_values = zeros(1, 59);
        sig_t_values(t_orig_t2(adj, task, :) > crit_t_t2(adj, task, 2)) = 1;
        sig_t_values(t_orig_t2(adj, task, :) < crit_t_t2(adj, task, 1)) = 1;
        sig_chan_number = find(sig_t_values == 1);

        for x = 1:length(sig_chan_number)
            sig_chans_all{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
        end

        figure; 
        if ~isempty(sig_chan_number)
            topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
        end
        topoplot(t_orig_t2(adj, task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
%         caxis([-7 0]); %c.Axis.FontSize = 16;
        colorbar;
    %     colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
        colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
        title([title_txt{task}, adj_text{adj}], 'FontSize', 30, 'FontWeight','normal')
        set(get(gca,'title'),'Position',[0,-.65, 0])
    end
    
    % testing with all participants together
    for task = 1:2
        [pval_t1_all(task, :), t_orig_t1_all(task, :), crit_t_t1_all(task, :),~,~]=mult_comp_perm_t1([coeffs{1, task}; coeffs{2, task}]);
    end

    % plot t-values from one sample t-test both groups together
    load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
    title_txt = {'Simple RT' 'Go/no-go'};%cmap = crameri('batlow');
    for task = 2%1:2

        sig_t_values = zeros(1, 59);
        sig_t_values(t_orig_t1_all(task, :) > crit_t_t1_all(task, 2)) = 1;
        sig_t_values(t_orig_t1_all(task, :) < crit_t_t1_all(task, 1)) = 1;
        sig_chan_number = find(sig_t_values == 1);

        for x = 1:length(sig_chan_number)
            sig_chans_all{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
        end

        figure; 
        if ~isempty(sig_chan_number)
            topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
        end
        topoplot(mean([coeffs{1, task}; coeffs{2, task}], 1), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        caxis([-.18 0]); %c.Axis.FontSize = 16;
        colorbar;
        colorbar('Ticks',[-.15, -.1, -.05, 0], 'FontSize', 26, 'FontWeight','normal');
        colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
%         title(title_txt{task}, 'FontSize', 30, 'FontWeight','normal')
%         set(get(gca,'title'),'Position',[0,-.65, 0])

    end
end

%% compare correlation coefficients with and without adjustment
% coeff_erp_pupil_phaseresid; % group x task
% coeff_erp_pupil;
clear sig_chans_all
for task = 1:2
    [pval_effect_of_adj(task, :), t_orig_effect_of_adj(task, :), crit_t_effect_of_adj(task, :),~,~] = mult_comp_perm_t1([coeff_erp_pupil_phaseresid{1, task}; coeff_erp_pupil_phaseresid{1, task}]-[coeff_erp_pupil{1, task}; coeff_erp_pupil{1, task}]);
    
    sig_t_values = zeros(1, 59);
    sig_t_values(t_orig_effect_of_adj(task, :) > crit_t_effect_of_adj(task, 2)) = 1;
    sig_t_values(t_orig_effect_of_adj(task, :) < crit_t_effect_of_adj(task, 1)) = 1;
    sig_chan_number = find(sig_t_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans_all{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end

    figure; 
    if ~isempty(sig_chan_number)
        topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
    end
    topoplot(t_orig_effect_of_adj(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([-4 0]); c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[-4, -2,  0], 'FontSize', 26, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
    title(title_txt{task}, 'FontSize', 30, 'FontWeight','normal')
    set(get(gca,'title'),'Position',[0,-.65, 0])
end

%% average correlation coefficients across channels and run repeated
% measures anova
load coeff_erp_pupil_phaseresid; % group x task
load coeff_erp_pupil;
for grp =1:2
    for task = 1:2
        coeffs_adj{grp, task} = mean(coeff_erp_pupil_phaseresid{grp, task}, 2);
        coeffs_not_adj{grp, task} = mean(coeff_erp_pupil{grp, task}, 2);
    end
end

% create table for spss
var = [[ones(length(coeffs_adj{1, 1}), 1); ones(length(coeffs_adj{2, 1}), 1)*2],...
    [coeffs_adj{1, 2}; coeffs_adj{2, 2}],[coeffs_not_adj{1, 2}; coeffs_not_adj{2, 2}]];

corr_r = array2table(var, "VariableNames",["group", "corr_r_adj" "cor_r_not_adj"]);

writetable(corr_r, 'corr_erp_pupil_phaseresid.xlsx')


%% correlation was strongest in which channel?
load coeff_erp_pupil_phaseresid
for grp = 1:2
    for task = 1:2
        [M(grp, task),I(grp, task)] = max(mean(coeff_erp_pupil_phaseresid{grp, task}, 1));
    end
end


%% plot C5 plot_all_data_2groups(data_grp1, data_grp2, y_label_text, title_text)
% coeff_robust_erpresid_chan {grp, task}(part, chan)
load coeff_robust_eegphaseresid; % correlation coefficients RT vs ERP residuals
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_txt = {'simpleRT', 'go/no-go'};
for c = 1:size(chanlocs_EEGChanOnly, 1)
    if strcmp(chanlocs_EEGChanOnly(c).labels, 'C5')
        chan = c;
    end
end
% plot_all_data_2tasks(data_grp1_task1, data_grp1_task2, data_grp2_task1, data_grp2_task2, y_label_text)
plot_all_data_2tasks(coeff_robust_eegphaseresid{1, 1}(:, chan), coeff_robust_eegphaseresid{1, 2}(:, chan),...
    coeff_robust_eegphaseresid{2, 1}(:, chan), coeff_robust_eegphaseresid{2, 2}(:, chan),'Correlation \itr')

% plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)
plot_all_data_onetask(coeff_robust_eegphaseresid{1, 2}(:, chan), coeff_robust_eegphaseresid{2, 2}(:, chan),'Correlation \itr')


%% stats using permutation method
load coeff_robust_eegphaseresid % correlation coefficients RT vs ERP residuals
clear pval_t1 t_orig_t1 crit_t_t1 pval_t1_all t_orig_t1_all crit_t_t1_all pval_t2 t_orig_t2 crit_t_t2
for task = 1:2
    for grp = 1:2
        [pval_t1(grp, task, :), t_orig_t1(grp, task, :), crit_t_t1(grp, task, :), est_alpha, seed_state]=mult_comp_perm_t1(coeff_robust_eegphaseresid{grp, task});
    end
    [pval_t2(task, :), t_orig_t2(task, :), crit_t_t2(task, :), est_alpha, seed_state] = mult_comp_perm_t2(coeff_robust_eegphaseresid{1, task},coeff_robust_eegphaseresid{2, task});
    [pval_t1_all(task, :), t_orig_t1_all(task, :), crit_t_t1_all(task, :), est_alpha, seed_state]=mult_comp_perm_t1([coeff_robust_eegphaseresid{1, task};coeff_robust_eegphaseresid{2, task}]);
end

% plot t-values from one sample t-test both groups together
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_txt = {'Simple RT' 'Go/no-go'};%cmap = crameri('batlow');
for task = 1:2
    
    sig_t_values = zeros(1, 59);
    sig_t_values(t_orig_t1_all(task, :) > crit_t_t1_all(task, 2)) = 1;
    sig_t_values(t_orig_t1_all(task, :) < crit_t_t1_all(task, 1)) = 1;
    sig_chan_number = find(sig_t_values == 1);

    for x = 1:length(sig_chan_number)
        sig_chans_all{task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
    end
    
%     figure; 
%     if ~isempty(sig_chan_number)
%         topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
%                 'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
%     end
%     topoplot(t_orig_t1_all(task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
%     caxis([0 5]); %c.Axis.FontSize = 16;
%     colorbar;
%     colorbar('Ticks',[0, 2, 4], 'FontSize', 30, 'FontWeight','normal');
%     colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
%     title(title_txt{task}, 'FontSize', 30, 'FontWeight','normal')
%     title('bamako', 'FontSize', 30, 'FontWeight','normal')
%     set(get(gca,'title'),'Position',[0,-.65, 0])


    % plot average correlation coefficients with significant electrodes
    % highlighted
    figure; 
    if ~isempty(sig_chan_number)
        topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
    end
    topoplot(mean([coeff_robust_eegphaseresid{1, task}; coeff_robust_eegphaseresid{2, task}], 1), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
    caxis([0 .1]); %c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[0, .05, .1], 'FontSize', 30, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php

end

%% plot electrodes where correlation r different from zero - t values
clear sig_chans
load G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs_EEGChanOnly.mat
title_txt = {'Young - simple RT' 'Older - simple RT'; 'Young - go/no-go' 'Older - go/no-go'};
for grp = 1:2
    for task = 1:2

        sig_t_values = zeros(1, 59);
        sig_t_values(t_orig_t1(grp, task, :) > crit_t_t1(grp, task, 2)) = 1;
        sig_t_values(t_orig_t1(grp, task, :) < crit_t_t1(grp, task, 1)) = 1;
        sig_chan_number = find(sig_t_values == 1);
       

        for x = 1:length(sig_chan_number)
            sig_chans{grp, task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
        end

        
        figure;
        if ~isempty(sig_chan_number)
            topoplot(ones(59, 1)*.025, chanlocs_EEGChanOnly, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
                'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
        end
        topoplot(t_orig_t1(grp, task, :), chanlocs_EEGChanOnly, 'electrodes', 'off'); 
        caxis([0 5]); %c.Axis.FontSize = 16;
        colorbar;
        colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
        colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
        title(title_txt{task, grp}, 'FontSize', 24, 'FontWeight','normal')
        set(get(gca,'title'),'Position',[0,-.6, 0])
%         text(5, 0.4, title_txt{task, grp})
        
    end
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
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    ax.Color = 'none';
%     ax.XTickLabel= [];
    ax.XTickLabel= [{'' 'Young' 'Older' ''}];
%     xticks([1 2])
    ylabel(y_label_text, 'FontSize', 32, 'FontWeight','normal')
    title(title_text, 'FontSize', 32, 'FontWeight','normal')
    
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
    xaxis = -.199:1/500:2;
    
    for x = 1.001:1/500:1.500 % grey background between 1000 amd 1500ms after ceu-onset
        plot([x x],[-12.8 1.8], 'color', [.8 .8 .8] ); hold on
    end
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
    axis([0 2 -13 2]);
    % title('GNG', 'FontSize', 24, 'FontWeight','bold')
    % legend('RT1', 'RT2', 'RT3', 'RT4', 'RT5', 'location', 'northwest')
    ax = gca;
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    title(title_text, 'FontSize', 32, 'FontWeight','normal');
    
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



function plot_all_data_onetask(data_grp1_task1, data_grp2_task1, y_label_text)

    % plot data for young group - go/nogo task
    figure; box on; hold on
    
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

