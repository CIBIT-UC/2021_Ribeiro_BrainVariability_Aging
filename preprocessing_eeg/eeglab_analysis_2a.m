% after determining ica weights for all experiment...
% Maria Ribeiro - May 2016
% manually select artifact components - only after epoching!
% epoch and baseline ica merged file and check ICs for artifact ICs and save their numbers for later removal
% save as artifact_ICs
% analyse individual runs separately for synchronization with eye data

clear
close all
participants=38;
task_load={'W1', 'D1', 'D2', 'G1', 'G2'};
task={'W1', 'D1', 'D2', 'G1', 'G2'};
% directory_ica='L:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ica_files';
directory_ica='N:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\ica_files';
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for p=participants
participant=strcat('AB', num2str(p));
directory=strcat('L:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\EEG\');

cd(directory)
filename=strcat('L:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\artifactICs\',participant, '_artifact_ICs.mat');
load(filename); load BadChannels; % for each participant

%clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
% load file with ica weights
if p==56
filename=strcat(participant, '_Merged_ica.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
else
filename=strcat(participant, '_Merged_ica_ChanLocs_epochs_baselined.set');
EEG = pop_loadset('filename',filename,'filepath',directory_ica);
end
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);

%% looping through the tasks
    for t=1:5
        % load individual file
        filename=strcat(participant, '_', task_load{t}, '_RefM1M2_RemEKG_EMG_Fps.set');
        EEG = pop_loadset('filename', filename,'filepath',directory);
        % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);

        % add measured channel locations
        filename=strcat('L:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\channel_locs\', participant, '_locs_corrected_61channels.DAT');
        EEG=pop_chanedit(EEG, 'load',filename);
        % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
        % filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_ChanLocs.set');
        % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',filename,'savenew',filename,'gui','off');

        %remove bad channels if necessary
        EEG = pop_select( EEG,'nochannel',BadChannels);
        % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
        % filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh.set');
        % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',filename,'savenew',filename,'gui','off'); 

        % add ica from dataset 1 into current dataset
        EEG = eeg_checkset( EEG );
        EEG = pop_editset(EEG, 'icachansind', 'ALLEEG(1).icachansind', 'icaweights', 'ALLEEG(1).icaweights', 'icasphere', 'ALLEEG(1).icasphere');
        % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
        % filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica.set');
        % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',filename,'savenew',filename,'gui','off'); 
        % 
        % filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica.set');
        % EEG = pop_loadset('filename', filename);

        % remove artifact components
        filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica_ArtICsRem');
        EEG = pop_subcomp( EEG, artifact_ICs, 0);
        % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
        % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',filename,'savenew',filename,'gui','off'); 

        % filter 0.1 - 35 Hz
        % EEG = eeg_checkset( EEG );
        filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica_ArtICsRem_Filt0_1_35Hz');
        EEG = pop_eegfiltnew(EEG, 0.1, 35, 16500, 0, [], 0); % filter order?
        % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',filename,'savenew',filename,'gui','off'); 
    end
end
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
eeglab redraw