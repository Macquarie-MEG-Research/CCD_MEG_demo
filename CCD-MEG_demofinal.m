
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script was created by Wei He and Robert Seymour for CCD-KIT MEG 
%                     Workshop on 6th Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Contents:

% - Step-0  : set-up
% - Step-1  : read the continuous data, segment into 2 seconds epochs
% - Step-2  : downsample the data to speed up component analysis
% - Step-3  : use ICA in order to identify cardiac and blink components
% - Step-4  : compute the power spectrum
% - Step-5  : compute the forward model
% - Step-6  : source reconstruction
% - Step-7  : create a pseudo-contrast based on a median split of the epoch
% - Step-8  : source reconstruction of the low-high alpha spectral data
% - Step-9  : connectivity analysis 
% - Step-10 : network analysis

%% Step-0: Set-Up

% You will first need to add Fieldtrip version v20171128 to your path. You
% only need to add the first folder - Fieldtrip will add the relevent bits
% of code from subfolders.

% This can be done using addpath('path_to_Fieldtrip')

% Now try typing ft_defaults into the Command Window. If there are no red
% error warnings you've done it correctly!

% % Below is the example script - just need to uncomment and change your
% % account name accordingly to make it work!
% if ispc 
%     addpath C:\Users\mq42604613\Desktop\fieldtrip-20171128
%     ft_defaults
%     clear all
% elseif ismac
%     addpath /Users/mq42604613/Desktop/fieldtrip-20171128
%     ft_defaults
%     clear all
% end
        
% You will also need to change your current directory (cd) to the location
% of the CCD_MEG_demo_final folder

% This can be done using cd('path_to_CD_MEG_demo_final')

%% Step-0.1: Visualising the MEG and head shape coregistration
load grad_trans.mat; % co-registered meg  
load sourcemodel.mat % load sourcemodel
load headmodel_singleshell.mat % load headmodel_singleshell


%%visualize the coregistration of sensors, and headshape.
figure;
ft_plot_sens(grad_trans, 'style', 'r*')
ft_plot_headshape(grad_trans.fid) %plot head 
view([0 -90 0])

%% Step-1: read the continuous data and segment into 2 seconds epochs
load data.mat; % individual meg data

cfg         = [];
cfg.length  = 2;
cfg.overlap = 0.5;
data        = ft_redefinetrial(cfg, data);

% this step is needed to 1) remove the DC-component, and to 
% 2) get rid of a few segments of data at the end of the recording, 
% which contains only 0's.

cfg        = [];
cfg.demean = 'yes';
cfg.trials = 1:(numel(data.trial)-6);
data       = ft_preprocessing(cfg, data);

%% Step-1.1:make a visual inspection and reject bad trials/sensors

% load yokogawa160.mat; % yokogawa 160 layout

cfg         = [];
cfg.method  = 'summary';
cfg.channel = 'MEG';
cfg.layout  = 'yokogawa160.mat'; 
dataclean   = ft_rejectvisual(cfg, data);

%%you can identify the rejected trial numbers with the followign code
badtrials = [];
for i=1:length(dataclean.cfg.artfctdef.summary.artifact)
  badtrials(i) = find(data.sampleinfo(:,1)==dataclean.cfg.artfctdef.summary.artifact(i));
end;
fprintf('\n Trial %d was rejected \n', badtrials);

%% Step-1.1:alternatively, you can use the list below
% This is the definition of the badtrials for the data that has been stored on disk:
badtrials  = [20	22	23	24	27	29	31	32	37	47	55	57	58	59	60	...
    68	74	89	93	94	121	122	137	138	143	154	156	162	167	168	170	173	174	...
    175	181	198	199	203	213	214	215	216	217	218	223	224	248	253	255	256	261	...
    266	275	283	286	287	288	289	291	302	303]; %KIT sample data 2718

cfg        = [];
cfg.trials = setdiff(1:numel(data.trial), badtrials);
dataclean  = ft_selectdata(cfg, data);

% STOP and double click the dataclean variable...
% You should now have 243 trials of data (2s long).

%% Step-2: downsample the data to speed up component analysis
dataclean.time(1:end) = dataclean.time(1); % this avoids numeric round off issues in the time axes upon resampling

cfg            = [];
cfg.resamplefs = 100; % Here we are downsampling from 1000Hz --> 100Hz
cfg.detrend    = 'yes'; % Helps with low-frequency drift
datads         = ft_resampledata(cfg, dataclean);

%% Setp-3: use ICA in order to identify cardiac and blink components
cfg                 = [];
cfg.method          = 'runica'; 
cfg.runica.maxsteps = 50; % here we are limiting the ICA decomposion to 50 steps for time
%cfg.randomseed     = 0;  % this can be uncommented to match the data that has been stored on disk
comp                = ft_componentanalysis(cfg, datads);

%%visualize components see: http://www.fieldtriptoolbox.org/example/use_independent_component_analysis_ica_to_remove_ecg_artifacts
cfg           = [];
cfg.component = [1:20];       % specify the component(s) that should be plotted
cfg.layout    = 'yokogawa160.mat'; % specify the layout file that should be used for plotting
cfg.viewmode = 'component';
% ft_topoplotIC(cfg, comp)
ft_databrowser(cfg, comp)

% these were the indices of the bad components that were identified 
% they may be different if you re-run the ICA decomposition with a random randomseed.
badcomp = [3 7 25]; % may be different numbers

cfg           = [];
cfg.component = badcomp;
dataica       = ft_rejectcomponent(cfg, comp);

% Alternatively you can decompose the original data as it was prior to 
% downsampling and project the ICs out of the data
% SEE:
% http://www.fieldtriptoolbox.org/example/use_independent_component_analysis_ica_to_remove_ecg_artifacts


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Exercise 1 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% finish the following exercise and try to see how data look different
% before and afetr preprocessing. 
%
% Using the ft_databrowser, you could also try to do an FFT to 
% compare the power spctrum pre- and post ICA

cfg          =[];
cfg.channel  = 'MEG';
cfg.viewmode = 'vertical';
figure;
ft_databrowser(cfg, data);
ft_databrowser(cfg, datads);
ft_databrowser(cfg, dataica);



%% Step-4.1: compute the power spectrum

% For more info see http://www.fieldtriptoolbox.org/tutorial/timefrequencyanalysis#time-frequency_analysis_iii

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'no';
datapow          = ft_freqanalysis(cfg, dataica);


%%plot the topography and the spectrum
figure;

cfg           = [];
cfg.layout    = 'yokogawa160.mat'; 
cfg.marker    = 'no'
cfg.xlim      = [9 11];
%cfg.colorbar  = 'East';
cfg.comment   = 'no'
subplot(211); ft_topoplotER(cfg, datapow);
title('Topoplot of 9-11Hz Power');
set(gca,'FontSize',10);

cfg         = [];
cfg.channel = {'AG130', 'AG131', 'AG132', 'AG134','AG137',...
    'AG145','AG146', 'AG148','AG150', 'AG153'};

subplot(212); ft_singleplotER(cfg, datapow);
xlabel('Frequency (Hz)');
title('Frequency Power Over 10 Occipital Sensors');
set(gca,'FontSize',10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Exercise 2 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following exercise investigates the effect of frequency smoothing on
% alpha power estimation

cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.pad    = 'nextpow2';
cfg.foilim = [0 50];
cfg.taper  = 'dpss';
cfg.tapsmofrq = 1;   % i.e. no real smoothing
freq       = ft_freqanalysis(cfg, dataica);
figure
semilogy(freq.freq, freq.powspctrm(1,:), 'b-');

cfg.taper     = 'dpss';
cfg.tapsmofrq = 3;   % i.e. lots of smoothing
freq          = ft_freqanalysis(cfg, dataica);
hold on
semilogy(freq.freq, freq.powspctrm(1,:), 'g-');
legend({'1Hz Smoothing', '3Hz Smoothing'});

xlabel('Frequency (Hz)');
ylabel('Power');
title({'The Effect of Smoothing with Multi-Tapers', 'on Frequency Analysis'});
set(gca,'FontSize',15);

%% Step-4.2: convert to planar and compute power spectrum

load neighbours.mat; % load  neighbourhood structure for the channels
 
% cfg                 = [];
% cfg.neighbours      = neighbours;
% cfg.layout    = 'yokogawa160.mat';
% ft_neighbourplot(cfg)

% convert axial to planar
dataicatmp      = dataica;
dataicatmp.grad = dataica.grad;

cfg               = [];
cfg.neighbours    = neighbours;
cfg.planarmethod  = 'sincos';
planar            = ft_megplanar(cfg, dataicatmp);

% compute the power spectrum 
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'no';
datapow_planar   = ft_freqanalysis(cfg, planar);

figure; 
cfg              = [];
cfg.layout       = 'yokogawa160.mat'; 
cfg.xlim         = [9 11];
subplot(2,2,1); ft_topoplotER(cfg, datapow); title('axial gradient representation')
subplot(2,2,2); ft_topoplotER(cfg, ft_combineplanar([], datapow_planar));title('planar gradient representation')

cfg         = [];
cfg.channel = {'AG130', 'AG131', 'AG132', 'AG134','AG137', 'AG145',...
    'AG146', 'AG148','AG150', 'AG153'};

subplot(2,2,[3 4]); ft_singleplotER(cfg, datapow);
title ('Power spectrum averaged over 10 occipital sensors ~10 Hz peak');
set(gca,'FontSize',10);


%% Step-5: compute the forward model

% compute the leadfield
cfg             = [];
cfg.grad        = grad_trans;
cfg.grid        = sourcemodel;
cfg.headmodel   = headmodel_singleshell;
cfg.reducerank  = 2; % default for MEG is 2, for EEG is 3
% cfg.normalize ='yes'; % control against the power bias towards the center of the head 
cfg.channel     = {'MEG'};
lf              = ft_prepare_leadfield(cfg,dataica);

% Make a figure to check everything is correct with the forward model
figure;
% make the headmodel surface transparent
ft_plot_vol(headmodel_singleshell,'facecolor', 'cortex', 'edgecolor', 'none'); camlight; alpha 0.4;           
ft_plot_mesh(lf.pos(lf.inside,:));
ft_plot_sens(grad_trans,'style','*r');
view([0 -90 0])

%% Step-6: source reconstruction
% compute sensor level Fourier spectra, to be used for cross-spectral density computation.
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;
cfg.pad        = 'nextpow2';
freq           = ft_freqanalysis(cfg, dataica);

% do the source reconstruction
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = lf;
cfg.headmodel         = headmodel_singleshell;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%'; % regularization parameter 
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.fixedori      = 'yes'; % only keep the largest of the three dipole directions per spatial filter 
sourceInt             = ft_sourcedescriptives([],ft_sourceanalysis(cfg, freq));

sourceInt.avg = rmfield(sourceInt.avg,{'ori','nai'}); % remove the orientation information due to FT bug

% Interpolate the solution onto an MRI overlay - you may experience some
% warnings - DON'T WORRY everything is working . Fieldtrip is just telling
% you that some (unused) fields are not the expected size. This is an annoying 
% bug and will probably change in future releases 

load mri_realigned.mat

cfg                  = [];
cfg.downsample       = 2;
cfg.parameter        = 'pow';
cfg.interpmethod     = 'nearest';
source_interpo       = ft_sourceinterpolate(cfg, sourceInt , mri_realigned);

% Plot this
cfg                  = [];
cfg.method           = 'ortho';
cfg.funparameter     = 'pow';
figure;
ft_sourceplot(cfg, source_interpo);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Exercise 3 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See how all the activity seems to be at the centre of the head! finish
% the exercise to find a simple solution. 

% Also see:
% http://www.fieldtriptoolbox.org/tutorial/beamformer#source_analysiswithout_contrasting_condition

%compute Neural-Activity-Index
sourceNAI = sourceInt;
% the power normalized with an estimate of the spatially inhomogeneous noise
sourceNAI.avg.pow = sourceInt.avg.pow ./ sourceInt.avg.noise;

cfg                  = [];
cfg.downsample       = 2;
cfg.parameter        = 'pow';
cfg.interpmethod = 'nearest';
sourceNAI_interpo       = ft_sourceinterpolate(cfg, sourceNAI , mri_realigned);

cfg                  = [];
cfg.method           = 'ortho';
cfg.funparameter     = 'pow';
figure;
ft_sourceplot(cfg, sourceNAI_interpo);

%% Step-7: creation of a pseudo-contrast based on a median split of the epoch
% compute sensor level single trial power spectra 
% NOTE: a similar step has been done in Step 4.1, what is the difference?

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [9 11];                          
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'yes';
cfg.pad          = 'nextpow2';
datapow          = ft_freqanalysis(cfg, dataica);

%%identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, 10);
tmp     = datapow.powspctrm(:,:,freqind); % trial (268)x channel (160) 
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max for the avraged trials
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));

%%compute the power spectrum for the median split data
cfg              = [];
cfg.trials       = indlow; 
datapow_low      = ft_freqdescriptives(cfg, datapow);

cfg.trials       = indhigh; 
datapow_high     = ft_freqdescriptives(cfg, datapow);

%%compute the difference between high and low alpha conditions
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'divide';
powratio      = ft_math(cfg, datapow_high, datapow_low);

%%plot the topography of the difference along with the spectra
figure;
cfg             = [];
cfg.layout      = 'yokogawa160.mat';
cfg.comment     = 'no';
cfg.xlim        = [9.9 10.1];
subplot(121); ft_topoplotTFR(cfg, powratio); 
title({'Difference topography of'; 'the median split data'})
set(gca,'FontSize',10);

%%compute the power spectrum in high and low-alpha trials
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';        
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'yes';
cfg.pad          = 'nextpow2';

cfg.trials       = indlow;      
datapow_low      = ft_freqanalysis(cfg, dataica);

cfg.trials       = indhigh;   
datapow_high      = ft_freqanalysis(cfg, dataica);

cfg         = [];
cfg.channel = powratio.grad.label(chanind);
subplot(122); ft_singleplotER(cfg, datapow_high, datapow_low);
xlabel('Frequency (Hz)');
legend({'High','Low'});
title({'Difference in power spectrum '; 'of the median split data'})
set(gca,'FontSize',10);


%%compute fourier spectra for frequency of interest according to the trial split
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;
cfg.pad        = 'nextpow2';

cfg.trials = indlow; 
freq_low   = ft_freqanalysis(cfg, dataica);

cfg.trials = indhigh; 
freq_high  = ft_freqanalysis(cfg, dataica);

%% Step-8: source reconstruction of the low-high alpha spectral data

%%compute the beamformer filters based on the entire data
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = lf;
cfg.headmodel         = headmodel_singleshell;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, freq);

% Why is it so important to use a 'common filter'?

% use the precomputed filters to localise power for low and high alpha epochs
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = lf;
cfg.grid.filter       = source.avg.filter;
cfg.headmodel         = headmodel_singleshell;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
source_low            = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_low));
source_high           = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_high));

source_low.avg = rmfield(source_low.avg,{'nai','csdlabel'}); % remove the CSD & NAI information due to FT bug
source_high.avg = rmfield(source_high.avg,{'nai','csdlabel'}); % remove the CSD & NAI information due to FT bug

% create a fancy mask to focus on the maximum differences in alpha power
cfg           = [];
cfg.operation = 'log10(x1)-log10(x2)';
cfg.parameter = 'pow';
source_ratio  = ft_math(cfg, source_high, source_low);
source_ratio  = rmfield(source_ratio,'freq'); % remove the freq information due to FT bug

source_ratio.mask    = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-1)))./2;

% Interpolate the data and mask onto the subject's MRI
cfg                  = [];
cfg.downsample       = 2;
cfg.parameter        = 'pow';
cfg.interpmethod     = 'nearest';
source_interpo       = ft_sourceinterpolate(cfg, source_ratio , mri_realigned);
cfg.parameter        = 'mask';
source_interpo_mask  = ft_sourceinterpolate(cfg, source_ratio, mri_realigned);
source_interpo.mask  = source_interpo_mask.mask;

% plot
cfg                  = [];
cfg.method           = 'ortho';
cfg.funparameter     = 'pow';
cfg.maskparameter    = 'mask';
cfg.location         = 'max';
ft_sourceplot(cfg, source_interpo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Exercise 4 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compare the source plot with sensor level plot and see how source map
% agree with speculations from the sensor plot 
%
% Can you also work out a fancier way to plot the results?
% http://www.fieldtriptoolbox.org/tutorial/plotting#plotting_data_at_the_source_level


%% Step-9: Connectivity analysis 
%%compute connectivity
cfg         = [];
cfg.method  ='coh';
cfg.complex = 'absimag'; % use the imaginary part of coherence
source_conn = ft_connectivityanalysis(cfg, source);
figure;imagesc(source_conn.cohspctrm(source_conn.inside, source_conn.inside));


%%source parcellation
load template_grid.mat; % template grid 
load atlas_mni.mat; % aal altas

source_conn.pos = template_grid.pos;

cfg                 = []; 
cfg.interpmethod    = 'nearest'; 
cfg.parameter       = 'tissue'; 
atlas_interp        = ft_sourceinterpolate(cfg, atlas, template_grid); 

atlas_interp.pos = template_grid.pos;

cfg                 = [];
cfg.parameter       = 'cohspctrm';
parc_conn           = ft_sourceparcellate(cfg, source_conn, atlas_interp);
figure;imagesc(parc_conn.cohspctrm);


%% Step-10: Network analysis

%full-connection
cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = .2;
network_fll  = ft_networkanalysis(cfg,source_conn);

network_full.pos = sourcemodel.pos;
network_full.dim = sourcemodel.dim;
network_full.inside = sourcemodel.inside;

%%plot
cfg                         = [];
cfg.downsample              = 2;
cfg.parameter               = 'degrees';
cfg.interpmethod            = 'nearest';
network_full_interpo        = ft_sourceinterpolate(cfg, network_full , mri_realigned);

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'degrees';
% cfg.funcolormap   = 'jet';
ft_sourceplot(cfg, network_full_interpo);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Exercise 5 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try a different threshold and see how the resutlt might differ?

% Try different connectivity esimates, eg., the phase locking value and
% see how the adjacency matrix and connectivity matrix may differ? 

