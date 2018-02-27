# mTRF_Toolbox
The mTRF Toolbox is a MATLAB repository that permits the fast computation of the linear stimulus-response mapping of any sensory system in the forward or backward direction. It is suitable for analysing multi-channel EEG, MEG, ECoG and EMG data. The forward encoding model, or temporal response function (TRF) as it is commonly known, can be used to investigate information processing in neuronal populations using conventional time-frequency and source analysis techniques. In addition, TRFs can be used to predict the spectrotemporal dynamics of future neural responses to unseen stimulus sequences. Stimulus reconstruction can also be performed using backward decoding models that project the multi-channel population responses back to the dynamics of the causal stimulus signal. The mTRF Toolbox facilitates the use of extended continuous stimuli in electrophysiological studies compared to conventional time-locked averaging approaches which require the use of discrete, isolated stimulus events. This allows researchers to investigate of how neural systems process dynamic environmental signals such as speech, music and motion.

Supporting documentation can be found in the accompanying open access paper:...
[1] Crosse MJ, Di Liberto GM, Bednar A, Lalor EC (2016) [The Multivariate Temporal Response Function (mTRF) Toolbox: A MATLAB Toolbox for Relating Neural Signals to Continuous Stimuli.](http://dx.doi.org/10.3389/fnhum.2016.00604) *Frontiers in Human Neuroscience* 10:604.

## Functions
- mTRFtrain.m
- mTRFpredict.m
- mTRFtransform.m
- mTRFcrossval.m
- mTRFmulticrossval.m
- lagGen.m

## Datasets
- speech_data.mat
- contrast_data.mat
- coherentmotion_data.mat
