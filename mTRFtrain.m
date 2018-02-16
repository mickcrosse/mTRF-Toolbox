function [model,t,c] = mTRFtrain(stim,resp,fs,map,tmin,tmax,lambda)
%mTRFtrain mTRF Toolbox training function.
%   MODEL = MTRFTRAIN(STIM,RESP,FS,MAP,TMIN,TMAX,LAMBDA) performs ridge
%   regression on the stimulus property STIM and the neural response data
%   RESP to solve for their linear mapping function MODEL. Pass in MAP==1
%   to map in the forward direction or MAP==-1 to map backwards. The
%   sampling frequency FS should be defined in Hertz and the time lags
%   should be set in milliseconds between TMIN and TMAX. Regularisation is
%   controlled by the ridge parameter LAMBDA.
%
%   [...,T,C] = MTRFTRAIN(...) also returns the vector of time lags T for
%   plotting MODEL and the regression constant C for absorbing any bias
%   when testing MODEL.
%
%   Inputs:
%   stim   - stimulus property (time by features)
%   resp   - neural response data (time by channels)
%   fs     - sampling frequency (Hz)
%   map    - mapping direction (forward==1, backward==-1)
%   tmin   - minimum time lag (ms)
%   tmax   - maximum time lag (ms)
%   lambda - ridge parameter
%
%   Outputs:
%   model  - linear mapping function (MAP==1: feats by lags by chans,
%            MAP==-1: chans by lags by feats)
%   t      - vector of time lags used (ms)
%   c      - regression constant
%
%   See README for examples of use.
%
%   See also LAGGEN MTRFTRANSFORM MTRFPREDICT MTRFCROSSVAL
%   MTRFMULTICROSSVAL.

%   References:
%      [1] Lalor EC, Pearlmutter BA, Reilly RB, McDarby G, Foxe JJ (2006)
%          The VESPA: a method for the rapid estimation of a visual evoked
%          potential. NeuroImage 32:1549-1561.
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2015) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Author: Edmund Lalor, Michael Crosse, Giovanni Di Liberto
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Email: edmundlalor@gmail.com
%   Website: www.lalorlab.net
%   April 2014; Last revision: Jan 8, 2016

% Define x and y
if tmin > tmax
    error('Value of TMIN must be < TMAX')
end
if map == 1
    x = stim;
    y = resp;
elseif map == -1
    x = resp;
    y = stim;
    [tmin,tmax] = deal(tmax,tmin);
else
    error('Value of MAP must be 1 (forward) or -1 (backward)')
end
clear stim resp

% Convert time lags to samples
tmin = floor(tmin/1e3*fs*map);
tmax = ceil(tmax/1e3*fs*map);

% Generate lag matrix
X = [ones(size(x)),lagGen(x,tmin:tmax)];

% Set up regularisation
dim = size(X,2);
if size(x,2) == 1
    d = 2*eye(dim,dim);d([1,end]) = 1;
    u = [zeros(dim,1),eye(dim,dim-1)];
    l = [zeros(1,dim);eye(dim-1,dim)];
    M = d-u-l;
else
    M = eye(dim,dim);
end

% Calculate model
model = (X'*X+lambda*M)\(X'*y);

% Format outputs
c = model(1:size(x,2),:);
model = reshape(model(size(x,2)+1:end,:),size(x,2),length(tmin:tmax),size(y,2));
t = (tmin:tmax)/fs*1e3;

end