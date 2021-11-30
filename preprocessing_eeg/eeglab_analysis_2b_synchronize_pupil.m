%% synchronize pupil data - task condition by task condition
clear; close all
participants=6;
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
task={'W1', 'D1', 'D2', 'G1', 'G2'};
task_load={'W1', 'D1', 'D2', 'G1', 'G2'};

for p=participants
participant=strcat('AB', num2str(p));
directory1=strcat('F:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\EEG\');
directory=strcat('F:\\ProjectAgingNeuromodulation\\AuditoryResearch\\AuditoryTask_EyeTracking_EEG_lab94\\', participant, '\\EEG\\');
cd(directory1)

%% clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
for t=1:5;
% load individual  eeg file
cd(directory1)
filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica_ArtICsRem_Filt0_1_35Hz.set');
EEG = pop_loadset('filename', filename,'filepath',directory1);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% import and synchronize pupil data - check which triggers are common to
% start and end synchronization - needs to be done manually so far!
filename_pupil=strcat('F:\\ProjectAgingNeuromodulation\\AuditoryResearch\\AuditoryTask_EyeTracking_EEG_lab94\\', participant,...
    '\\ET\\', participant, '_', task_load{t},'_BlinkCorrectedWithEvents.mat');
if t==1;
    EEG = pop_importeyetracker(EEG,filename_pupil,[1 10] ,[5:10] ,{'R_Dia_X_(mm)' 'R_Dia_Y_(mm)' 'R_Dia_X_BCS_(mm)' 'R_Dia_Y_BCS_(mm)' 'R_Dia_X_BCU_(mm)' 'R_Dia_Y_BCU_(mm)'},0,1,0,1);
else
EEG = pop_importeyetracker(EEG,filename_pupil,[9 10] ,[5:10] ,{'R_Dia_X_(mm)' 'R_Dia_Y_(mm)' 'R_Dia_X_BCS_(mm)' 'R_Dia_Y_BCS_(mm)' 'R_Dia_X_BCU_(mm)' 'R_Dia_Y_BCU_(mm)'},0,1,0,1);
end
% save data
filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_ChanLocs_RemBadCh_ica_ArtICsRem_Filt0_1_35Hz_Pupil');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'savenew',filename,'gui','off'); 
% EEG = eeg_checkset( EEG );
% pop_eegplot( EEG);
end
end
eeglab redraw
