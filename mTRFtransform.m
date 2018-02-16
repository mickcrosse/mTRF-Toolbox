function [model_t,t,c_t] = mTRFtransform(stim,resp,model,fs,map,tmin,tmax,c)
%mTRFtransform mTRF Toolbox mapping transformation function.
%   MODEL_T = MTRFTRANSFORM(B_MODEL,RESP,STIM) tansforms the coefficients
%   of the model weights MODEL into transformed model coefficients MODEL_T.
%
%   Inputs:
%   stim   - stimulus property (time by features)
%   resp   - neural response data (time by channels)
%   model  - linear mapping function (MAP==1: feats by lags by chans,
%            MAP==-1: chans by lags by feats)
%   fs     - sampling frequency (Hz)
%   map    - transformation direction (forward -> backward==1, backward ->
%            forward==-1)
%   tmin   - minimum time lag (ms)
%   tmax   - maximum time lag (ms)
%   c      - regression constant
%
%   Outputs:
%   model_t - transformed model weights (lags by chans)
%   t       - vector of time lags used (ms)
%   c_t     - transformed model constant
%
%   See README for examples of use.
%
%   See also LAGGEN MTRFTRAIN MTRFPREDICT MTRFCROSSVAL MTRFMULTICROSSVAL.

%   References:
%      [1] Haufe S, Meinecke F, Gorgen K, Dahne S, Haynes JD, Blankertz B,
%          Bießmann F (2014) On the interpretation of weight vectors of
%          linear models in multivariate neuroimaging. NeuroImage 87:96-110.

%   Authors: Adam Bednar, Emily Teoh, Giovanni Di Liberto, Michael Crosse
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Email: edmundlalor@gmail.com
%   Website: www.lalorlab.net
%   April 2016; Last revision: 15 July 2016

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

% Convert time lags to samples
tmin = floor(tmin/1e3*fs*map);
tmax = ceil(tmax/1e3*fs*map);

% Generate lag matrix
X = [ones(size(x)),lagGen(x,tmin:tmax)];

% Transform model weights
model = [c;reshape(model,size(model,1)*size(model,2),size(model,3))];
model_t = (X'*X)*model*inv(y'*y);

% Format outputs
c_t = model_t(1:size(x,2),:);
model_t = reshape(model_t(size(x,2)+1:end,:),size(x,2),length(tmin:tmax),size(y,2));
t = (tmin:tmax)/fs*1e3;

end