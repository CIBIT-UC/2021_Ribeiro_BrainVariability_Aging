% interpolating blinks Maria Ribeiro 2016
clear
participant='AB9';
task={'W1', 'D1', 'D2', 'G1', 'G2'};
% directory=strcat('E:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\', participant, '\ET\');
directory=pwd;
cd(directory);

for t=2%1:5;
clearvars -except t participant task directory
filename=strcat(participant, '_', task(t), 'Events');
load(filename{1})
pupil_X=data(:,5)';
pupil_Y=data(:,6)';
graphics=1;
lrtask=0;
% manualblinks
% lowthresh
BlinkCorrection_X=stublinks(pupil_X,graphics,lrtask);
title_text=strcat(task(t), 'pupil X smoothed');
figure; plot(BlinkCorrection_X.NoBlinks); axis([-inf inf 1 7]);title(title_text{1});
title_text=strcat(task(t),'pupil X unsmoothed');
figure; plot(BlinkCorrection_X.NoBlinksUnsmoothed); axis([-inf inf 1 7]); title(title_text{1})
% BlinkCorrection=ans;
% figure; plot(BlinkCorrection.NoBlinksUnsmoothed)
BlinkCorrection_Y=stublinks(pupil_Y,graphics,lrtask);  
title_text=strcat(task(t),'pupil Y smoothed');
figure; plot(BlinkCorrection_Y.NoBlinks); axis([-inf inf 1 7]); title(title_text{1})
title_text=strcat(task(t),'pupil Y unsmoothed');
figure; plot(BlinkCorrection_Y.NoBlinksUnsmoothed); axis([-inf inf 1 7]); title(title_text{1})

data(:, 7)=BlinkCorrection_X.NoBlinks';
data(:, 8)=BlinkCorrection_Y.NoBlinks';
colheader{7}='R Dia X BCS [mm]';
colheader{8}='R Dia Y BCS [mm]';

data(:, 9)=BlinkCorrection_X.NoBlinksUnsmoothed';
data(:, 10)=BlinkCorrection_Y.NoBlinksUnsmoothed';
colheader{9}='R Dia X BCU [mm]';
colheader{10}='R Dia Y BCU [mm]';

filename_save=strcat(participant, '_', task(t), '_BlinkCorrectedWithEvents.mat');
% save(filename_save{1}, 'colheader', 'comments', 'data', 'messages', 'event');

%% saving data file to import into eeglab
PupilData=data(:,[5:10 12]);
filename_save=strcat(participant, '_', task(t), '_PupilData.mat');
% save(filename_save{1}, 'PupilData');
end






