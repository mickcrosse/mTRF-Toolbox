function [pred,r,p,rmse] = mTRFpredict(stim,resp,model,fs,map,tmin,tmax,c)
%mTRFpredict mTRF Toolbox prediction function.
%   PRED = MTRFPREDICT(STIM,RESP,MODEL,FS,MAP,TMIN,TMAX,C) performs a
%   convolution of the stimulus property STIM or the neural response data
%   RESP with their linear mapping function MODEL to solve for the
%   prediction PRED. Pass in MAP==1 to predict RESP or MAP==-1 to predict
%   STIM. The sampling frequency FS should be defined in Hertz and the time
%   lags should be set in milliseconds between TMIN and TMAX. The
%   regression constant C absorbs any bias in MODEL.
%
%   [...,R,P,RMSE] = MTRFPREDICT(...) also returns the correlation
%   coefficients R between the original and predicted values, the
%   corresponding p-values P and the root-mean-square errors RMSE.
%
%   Inputs:
%   stim   - stimulus property (time by features)
%   resp   - neural response data (time by channels)
%   model  - linear mapping function (MAP==1: feats by lags by chans,
%            MAP==-1: chans by lags by feats)
%   fs     - sampling frequency (Hz)
%   map    - mapping direction (forward==1, backward==-1)
%   tmin   - minimum time lag (ms)
%   tmax   - maximum time lag (ms)
%   c      - regression constant
%
%   Outputs:
%   pred   - prediction (MAP==1: time by chans, MAP==-1: time by feats)
%   r      - correlation coefficients
%   p      - p-values of the correlations
%   rmse   - root-mean-square errors
%
%   See README for examples of use.
%
%   See also LAGGEN MTRFTRAIN MTRFCROSSVAL MTRFMULTICROSSVAL.

%   Authors: Mick Crosse, Giovanni Di Liberto, Edmund Lalor
%   Email: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Website: www.lalorlab.net
%   Lalor Lab, Trinity College Dublin, IRELAND
%   April 2014; Last revision: 4-Feb-2019

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
X = [ones(size(x,1),1),lagGen(x,tmin:tmax)];

% Calculate prediction
model = [c;reshape(model,size(model,1)*size(model,2),size(model,3))];
pred = X*model;

% Calculate accuracy
if ~isempty(y)
    [r,p] = corr(y,pred);
    r = diag(r); p = diag(p);
    rmse = sqrt(mean((y-pred).^2,1))';
end

end