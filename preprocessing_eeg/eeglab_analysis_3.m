% STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
clear; close all;
participants=6;
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for p=participants;
participant=strcat('AB', num2str(p));
task={'W1', 'D1', 'D2', 'G1', 'G2'};
directory1=strcat('F:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\EEG\');
directory=strcat('F:\\ProjectAgingNeuromodulation\\AuditoryResearch\\AuditoryTask_EyeTracking_EEG_lab94\\', participant, '\\EEG\\');

cd(directory1);
% data files from eeglab_analysis_2b_synchronize_pupil
% _RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica_ArtICsRem_Filt0_1_35Hz_Pupil.set
% manual/visual removal of periods with artifacts
% save as ..._Filt0_1_35Hz_Pupil_ManualArtRej.set

% epoch data and baseline
for t=1:5
% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

% detection task
for t=2:3;
%load file
filename=strcat(participant, '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej.set');
EEG = pop_loadset(filename);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% epoch data
filename=strcat(participant, '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej_epochs');
EEG = pop_epoch( EEG, {  '1'  }, [-0.2           6], 'newname', filename, 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew',filename,'gui','off'); 
EEG = eeg_checkset( EEG );
% baseline epochs
EEG = pop_rmbase( EEG, [-200    0]);
filename=strcat(participant, '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej_epochs_baselined');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',filename,'savenew',filename,'gui','off');
end

% merge task runs
% name it as _Merged_epochs_baselined
% detection task
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
filename=strcat(participant, '_', task{2}, '_Filt0_1_35Hz_Pupil_ManualArtRej_epochs_baselined.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
filename=strcat(participant, '_', task{3}, '_Filt0_1_35Hz_Pupil_ManualArtRej_epochs_baselined.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
EEG = pop_mergeset( ALLEEG, [1  2], 1);
filename=strcat(participant, '_D_Merged_epochs_baselined');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',filename,'savenew',filename,'gui','off'); 
%% gng task
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
filename=strcat(participant, '_', task{4}, '_Filt0_1_35Hz_Pupil_ManualArtRej_epochs_baselined.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
filename=strcat(participant, '_', task{5}, '_Filt0_1_35Hz_Pupil_ManualArtRej_epochs_baselined.set');
EEG = pop_loadset('filename',filename,'filepath',directory);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
EEG = pop_mergeset( ALLEEG, [1  2], 1);
filename=strcat(participant, '_G_Merged_epochs_baselined');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',filename,'savenew',filename,'gui','off'); 

end
eeglab redraw