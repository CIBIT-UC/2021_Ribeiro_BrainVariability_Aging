% calculating behavior variables per run in epoched data before artifact rejection
% 15 March 2017 - Maria Ribeiro
clear
younger=[4	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46 50 51 53 54 56	59	62 66	68	72	74	76  78  80 81   82  84  85];
older=[7  8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77 79   83  86];

participants=[younger older];
task={'W1', 'D1', 'D2', 'G1', 'G2'};

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for p=participants
%clear eeglab
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    clearvars -except p b participants task test_errors2 ALLEEG EEG CURRENTSET ALLCOM LASTCOM PLUGINLIST STUDY
%  close all
directory=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\');
cd(directory);
load rejected_epochs;

    % detection and gng task
    for t=2:5

        % load behavior and calculate error trials and trials after error 
        cd(strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\',  task{t},'_behavior'));
        clear misses multiple_responses response2cue slow_responses response2nogo error_trials RejectEpochs;
        % error_trials - all errors, response2cue, multiple responses, misses, slow responses
        load correct_trial.mat;
        load misses.mat;
        load multiple_responses.mat;
        load response2cue.mat;
        load slow_responses.mat;

        if t==2 || t==3
            error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses);
        else 
            load response2nogo.mat
            error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses, response2nogo);    
        end
        
        % variable with epochs to reject including artifact epochs and
        % errors and trials after errors
        RejectEpochsCorrect{t-1}=sortrows(unique(cat(1, setdiff(1:60,correct_trial)', rejected_epochs{t}', [error_trials; error_trials+1])));
        % in the case of the last trial being an error need to delete the error trial +1 as it would be outside the existign trials
        RejectEpochsCorrect{t-1}(RejectEpochsCorrect{t-1}==61)=[];
            
    % back to ET directory
    cd(directory);
        % clear eeglab
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        % load file after baseline correction
        filename=strcat('AB', num2str(p),  '_', task(t), '_epochs_baselined.set');
        EEG = pop_loadset('filename',filename,'filepath',directory);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        % remove bad trials, error trials, and trials after errors and cue-only trials and no-go trials!
        EEG = pop_select( EEG,'notrial',RejectEpochsCorrect{t-1});

        % save variable with pupil epoched data X smoothed data
        PupilDilation=squeeze(EEG.data(3, :,:));     PupilDilation=PupilDilation';

        % load file before baseline correction
        filename=strcat('AB', num2str(p),  '_', task(t), '_epochs.set');
        EEG = pop_loadset('filename',filename,'filepath',directory);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        % remove bad trials, error trials, and trials after errors and cue-only trials and no-go trials!
        EEG = pop_select( EEG,'notrial',RejectEpochsCorrect{t-1});

        % save variable with pupil epoched data X smoothed data
        PupilDilation_Raw=squeeze(EEG.data(3, :,:));     PupilDilation_Raw=PupilDilation_Raw';
        PupilBaseline=mean(PupilDilation_Raw(:,192:240), 2);

        % create variable with cue-locked pupil response divided by pupil baseline of that trial, trial-by-trial
        % as pupil data has already been subtracted by baseline this further
        % division by baseline creates variable with percentage of change from
        % baseline PercChange=(PR-PB)/PB*100
        % where PB is pupil baseline before baseline correction, obviously!
        clear BaselineDivided_PupilDilation_Correct;
            for i=1: size(PupilDilation, 1);
                BaselineDivided_PupilDilation_Correct(i, :)=PupilDilation(i, :)/PupilBaseline(i)*100;
            end

        BslnDivided_PupilDilation_CorrectPerRun{t-1, :, :}=BaselineDivided_PupilDilation_Correct;
        
        % include reaction time variable for correlations
        clear epoch_event;
epoch_event=zeros(size(EEG.event, 2), 3);

for e=1:size(EEG.event, 2);
    epoch_event(e, 1)=EEG.event(e).epoch;
    epoch_event(e, 2)=EEG.event(e).type;
    epoch_event(e, 3)=EEG.event(e).latency;
end

detection_RT=[];
warning2imperativetime_detection=[];

        for r=1:epoch_event(end, 1);
            events_per_epoch_index=[];
            events_per_epoch_index=find(epoch_event(:,1)==r);

            clear events_per_epoch latencies_per_epoch
            events_per_epoch=epoch_event(events_per_epoch_index, 2);
            latencies_per_epoch=epoch_event(events_per_epoch_index, 3);
    
            if events_per_epoch(2)==2 && events_per_epoch(3)==5 && length(events_per_epoch)==3
                detection_RT=[detection_RT; r (latencies_per_epoch(3)-latencies_per_epoch(2))*1000/240];
                warning2imperativetime_detection=[warning2imperativetime_detection; r (latencies_per_epoch(2)-latencies_per_epoch(1))*1000/240]; % added 19 September 2016
            end  
        end 

        ReactionTime_CorrectPerRun{t-1, :, :}=[detection_RT warning2imperativetime_detection(:,2)];
    end
save BslnDivided_PupilDilation_CorrectPerRun BslnDivided_PupilDilation_CorrectPerRun;
save RejectEpochsCorrect RejectEpochsCorrect
save ReactionTime_CorrectPerRun ReactionTime_CorrectPerRun;
end
