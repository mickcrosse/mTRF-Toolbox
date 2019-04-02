# mTRF-Toolbox

The mTRF Toolbox is a MATLAB repository that permits the fast computation 
of the linear stimulus-response mapping of any sensory system in the 
forward or backward direction. It is suitable for analysing multi-channel 
EEG, MEG, ECoG and EMG data. The forward encoding model, or temporal 
response function (TRF) as it is commonly known, can be used to investigate 
information processing in neuronal populations using conventional time-
frequency and source analysis techniques. In addition, TRFs can be used to 
predict the spectrotemporal dynamics of future neural responses to unseen 
stimulus sequences. Stimulus reconstruction can also be performed using 
backward decoding models that project the multi-channel population 
responses back to the dynamics of the causal stimulus signal. The mTRF 
Toolbox facilitates the use of extended continuous stimuli in 
electrophysiological studies compared to conventional time-locked averaging 
approaches which require the use of discrete, isolated stimulus events. 
This allows researchers to investigate of how neural systems process 
dynamic environmental signals such as speech, music and motion.

## Tips on Practical Use

* Ensure that the stimulus and response data have the same sample rate
  and number of samples.
* Downsample the data when conducting large-scale multivariate analyses
  to reduce running time, e.g., 128 Hz or 64 Hz.
* Normalise or standardise the data beforehand. We recommend normalising 
  by the standard deviation. This will stabalise regularization across 
  trials/subjects/groups and facilitate a smaller parameter search.
* Enter the start and finish time lags in milliseconds. Enter positive
  lags for post-stimulus mapping and negative lags for pre-stimulus
  mapping. This is the same for both forward and backward mapping - the 
  code will automatically reverse the lags for backward mapping.
* When using MTRFPREDICT, always enter the model in its original
  3-dimensional form, i.e., do not remove any singleton dimensions.
* When using MTRFCROSSVAL, the trials do not have to be the same length,
  but using trials of the same length will optimise performance.
* When using MTRFMULTICROSSVAL, the trials in each of the three sensory
  conditions should correspond to the stimuli in STIM.

## Additional Information

mTRF Toolbox is also available for download at:
SourceForge: http://sourceforge.net/projects/aespa

mTRF Toolbox support documentation is available at:
http://dx.doi.org/10.3389/fnhum.2016.00604

For any questions and comments, please email:
mickcrosse@gmail.com (Mick Crosse) or edmundlalor@gmail.com (Ed Lalor)
