function [model,t] = mTRFtrain(stim,resp,fs,map,tmin,tmax,lambda,tlims)
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
%   stim   - stimulus property (time by features), this can be a cell array
%            where each cell is a different trial
%   resp   - neural response data (time by channels), this can be a cell array
%            where each cell is a different trial
%   fs     - sampling frequency (Hz)
%   map    - mapping direction (forward==1, backward==-1)
%   tmin   - minimum time lag (ms)
%   tmax   - maximum time lag (ms)
%   lambda - ridge parameter
%   tlims  - (optional) (NEW, NZ, 2019) specifies range or indexes of times 
%      that should be included in the model training and testing. If specific 
%      indexes are desired, then they should be specified in each cell of 
%      tlims, where the number of cells equals the number of trials.
%      (default: all indexes are used)
%      (see usetinds.m for more information on specifying tlims)
%
%   Outputs:
%   model  - linear mapping function (MAP==1: feats by lags by chans,
%            MAP==-1: chans by lags by feats)
%   t      - vector of time lags used (ms)
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

%   Authors: Mick Crosse, Giovanni Di Liberto, Edmund Lalor
%   Email: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Website: www.lalorlab.net
%   Lalor Lab, Trinity College Dublin, IRELAND
%   April 2014; Last revision: 4-Feb-2019

% If tlims isn't specified, use all indexes
if nargin<8, tlims = []; end

% If x and y are not cell arrays (if they contain data from just one trial,
% for example), make them cell arrays with one cell
% (needed to make the design matrices later)
if ~iscell(x), x = {x}; end
if ~iscell(y), y = {y}; end

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

% Generate lag matrix for each trial
%X = [ones(size(x,1),1),lagGen(x,tmin:tmax)];
disp('Creating the design matrices...');
% Create the design matrices trial by trial
std_tm = tic;
for i = 1:numel(x)
    % Generate lag matrix
    x{i} = [ones(size(x{i},1),1),lagGen(x{i},tmin:tmax)]; %%% skip constant term (NZ)
    % Set X and y to the same length
    minlen = min([size(x{i},1) size(y{i},1)]);
    x{i} = x{i}(1:minlen,:);
    y{i} = y{i}(1:minlen,:);
    % Remove time indexes, if specified
    if iscell(tlims), % if tlims is a cell array, it means that specific indexes were supplied
        tinds = tlims{i};
    else
        tinds = usetinds(tlims,fs,minlen);
    end
    x{i} = x{i}(tinds,:);
    y{i} = y{i}(tinds,:);
end
fprintf('Completed in %.3f s\n',toc(std_tm));

% Set up regularisation
dim = size(X,2);
if size(x,2) == 1
    d = 2*eye(dim);d([dim+2,end]) = 1;
    u = [zeros(dim,1),eye(dim,dim-1)];
    l = [zeros(1,dim);eye(dim-1,dim)];
    M = d-u-l; M(:,1) = 0; M(1,:) = 0;
else
    M = eye(dim,dim); M(1,1) = 0;
end

fprintf('Training the model...');
mdl_tm = tic;
[xtx,xty] = compute_linreg_matrices(x,y);
% Calculate model
% model = (X'*X+lambda*M)\(X'*y);
model = (xtx+lambda*M)\xty;
fprintf('Completed in %.3f s\n',toc(mdl_tm));

% Format outputs
t = (tmin:tmax)/fs*1e3;
% c = model(1,:);
% model = reshape(model(2:end,:),size(x,2),length(t),size(y,2));
model = reshape(model,size(x,2),length(t),size(y,2));

end