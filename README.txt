Code for

Slow fluctuations in ongoing brain activity decrease in amplitude with ageing yet their impact on task-related evoked responses is dissociable from behaviour
Maria J. Ribeiro, Miguel Castelo-Branco
bioRxiv 2021.11.18.467545; doi: https://doi.org/10.1101/2021.11.18.467545 


Task and acquisition parameters are described in detail in
Age-related differences in event-related potentials and pupillary responses in cued reaction time tasks.
Maria J. Ribeiro, Miguel Castelo-Branco. Neurobiol Aging. 2019 Jan;73:177-189. doi: 10.1016/j.neurobiolaging.2018.09.028


Analyses scripts run on Matlab and depend on:

	eeglab functions - https://sccn.ucsd.edu/eeglab/index.php
		in topoplot function change line 275 to COLORARRAY  = { [0 0 0] [0.5 0 0] [0 0 0] }; to make circles black
	
	Pupillary waveform deblinking Matlab code stublinks.m from https://sites.pitt.edu/~gsiegle/

	colour maps from http://www.fabiocrameri.ch/colourmaps.php
	https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

	permutation stats from David Groppe:
	https://www.mathworks.com/matlabcentral/fileexchange/29782-mult_comp_perm_t1-data-n_perm-tail-alpha_level-mu-reports-seed_state
	https://www.mathworks.com/matlabcentral/fileexchange/54585-mult_comp_perm_t2-data1-data2-n_perm-tail-alpha_level-mu-t_stat-reports-seed_state

	mintersect - https://www.mathworks.com/matlabcentral/fileexchange/6144-mintersect-multiple-set-intersection

	jbfill - https://www.mathworks.com/matlabcentral/fileexchange/13188-shade-area-between-two-curves

	fooof toolbox – https://fooof-tools.github.io/fooof/index.html

	Robust Correlation Toolbox v2 - https://github.com/CPernet/Robust-Correlations


Task scripts depend on:

	Psychophysics Toolbox Version 3 - http://psychtoolbox.org/


Raw EEG and pupil data available in OpenNeuro Repository - https://openneuro.org/datasets/ds003690/versions/1.0.0
Maria J. Ribeiro and Miguel Castelo-Branco (2021). EEG, ECG and pupil data from young and older adults: rest and auditory cued reaction time tasks. OpenNeuro. 
[Dataset] doi: 10.18112/openneuro.ds003690.v1.0.0


For questions, feel free to contact by email mjribeiro@fmed.uc.pt; twitter @Ribeiro_neurosc