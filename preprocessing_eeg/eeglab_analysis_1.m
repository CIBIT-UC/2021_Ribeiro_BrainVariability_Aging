% data analysis for running ica on merged files
% Maria Ribeiro May 2016
% start eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];

participant='AB6';
task={'W1', 'D1', 'D2', 'G1', 'G2'};
% task2=[1 3 4 5];
% for p=1;
    for t=1:5
% import .cnt file
filename=fullfile(pwd, strcat(participant,'_', task{t}, '.cnt'));
EEG = pop_loadcnt(filename);
% filename=strcat(participant, '_', task{t});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',filename,'gui','off'); 
EEG = eeg_checkset( EEG );
% save .set file
filename=strcat(participant, '_', task{t}, '.set');
EEG = pop_saveset( EEG, 'filename',filename,'filepath',pwd);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% re-reference to linked M1 and M2 earlobes
EEG = pop_reref( EEG, [33 43] ,'exclude',[65 66 67 68] );
filename=strcat(participant, '_', task{t}, '_RefM1M2');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'savenew',filename,'gui','off'); 

% remove EKG and EMG and Fps
EEG = pop_select( EEG,'nochannel',{'FP1' 'FPZ' 'FP2' 'EKG' 'EMG'});
filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',filename,'savenew',filename,'gui','off'); 

% filter eeg data between 1 and 100 Hz
EEG = pop_eegfiltnew(EEG, 1, 100, 1650, 0, [], 0);
filename=strcat(participant, '_', task{t}, '_RefM1M2_RemEKG_EMG_Fps_filtered1_100Hz');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'savenew',filename,'gui','off'); 
    end
% end

eeglab redraw

% manually reject bad epochs and bad channels - save as _ManualArtRej
% run ica on merged files

