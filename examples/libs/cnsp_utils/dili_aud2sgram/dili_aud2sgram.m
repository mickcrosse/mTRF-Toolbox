function [env,sgram,f,t] = dili_aud2sgram(filename,downFs)
% This function returns the Hilbert envelope and spectrogram (including the
% time and frequency axes) from a given file at the downsampling frequency
% 'downFs'.
% Spectrograms are extracted in 8 frequency bands according to the
% Greenwood equation. 'f' is the start frequency for each bin and 't' the
% time vector in seconds.
% ** Note that only the first audio channel is used in this function.
%
% Example: [env,sgram,f,t] = dili_aud2sgram('audio1.wav',100);
%
% This function uses the myspectrogram function (Kamil Wojcicki)
%
% Di Liberto-lab. Last update: 9 March 2023
%

% Eight frequency bands according to the Greenwood equation
greenwoodFreq = [250,     437.3751832
                437.3751832, 709.2700085;...
                709.2700085, 1103.808894;...
                1103.808894, 1676.313105;...
                1676.313105, 2507.057762;...
                2507.057762, 3712.527821;...
                3712.527821, 5461.751307;...
                5461.751307, 8000];

% Loading audio
[audio,fs] = audioread(filename);
audio = audio(:,1);

% Calculating envelope
env = abs(hilbert(audio));
if mod(fs,downFs) == 0 % if integer downsampling factor
    env = downsample(env,fs/downFs);
else
    env = resample(env,downFs,fs,0);
end

% Getting spectrogram
[~,sgramF,sgramData] = myspectrogram(audio,fs,[18 1], @hamming, 1024, [-45 -2], [1 -0.97], 'default', false, 'per');

% Make 8 frequency bands (logarithmically scaled), as defined by the
% Greenwood Equation, out of the original 1024 frequency bands
for ffreq = 1:size(greenwoodFreq, 1)
    freqIdxs = sgramF >= greenwoodFreq(ffreq,1) & sgramF < greenwoodFreq(ffreq,2);
    sgram(:,ffreq) = mean(sgramData(freqIdxs,:),1)';
end

% Downsampling sgram
for fband = 1:size(sgram,2)
    sgramDown(:,fband) = resamplee(sgram(:,fband),length(env),size(sgram,1));
end
sgram = sgramDown;
f = greenwoodFreq(:,1);
durSec = length(env)/downFs;
t = 0:(durSec/(length(env)-1)):durSec;

end

% Dealing with difficult downsampling factors
function resig = resamplee(sig,upsample,downsample)
    if upsample*downsample<2^31
        resig = resample(sig,upsample,downsample,0); % no filtering applied
    else
        sig1half = sig(1:floor(length(sig)/2));
        sig2half = sig(floor(length(sig)/2):end);
        resig1half = resamplee(sig1half,floor(upsample/2),length(sig1half));
        resig2half = resamplee(sig2half,upsample-floor(upsample/2),length(sig2half));
        resig = [resig1half;resig2half];
    end
end
