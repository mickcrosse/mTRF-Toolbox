function [pred,r,p,rmse] = mTRFpredict(stim,resp,model,fs,map,tmin,tmax,tlims)
%mTRFpredict mTRF Toolbox prediction function.
%   PRED = MTRFPREDICT(STIM,RESP,MODEL,FS,MAP,TMIN,TMAX) performs a
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
%   stim   - stimulus property (time by features), this can be a cell array
%            where each cell is a different trial
%   resp   - neural response data (time by channels), this can be a cell array
%            where each cell is a different trial
%   model  - linear mapping function (MAP==1: feats by lags by chans,
%            MAP==-1: chans by lags by feats)
%   fs     - sampling frequency (Hz)
%   map    - mapping direction (forward==1, backward==-1)
%   tmin   - minimum time lag (ms)
%   tmax   - maximum time lag (ms)
%   tlims  - (optional) (NEW, NZ, 2019) specifies range or indexes of times 
%      that should be included in the model training and testing. If specific 
%      indexes are desired, then they should be specified in each cell of 
%      tlims, where the number of cells equals the number of trials.
%      (default: all indexes are used)
%      (see usetinds.m for more information on specifying tlims)
%
%   Outputs:
%   pred   - prediction (MAP==1: time by chans, MAP==-1: time by feats)
%   r      - correlation coefficients
%   p      - p-values of the correlations
%   rmse   - root-mean-square errors
%     ** For all measures of prediction accuracy (r, p, mse):
%          MAP==1: trial by feats
%          MAP==-1: trial by chans
%
%   See README for examples of use.
%
%   See also LAGGEN MTRFTRAIN MTRFCROSSVAL MTRFMULTICROSSVAL.

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
if ~iscell(stim), stim = {stim}; end
if ~iscell(resp), resp = {resp}; end

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
% X = [ones(size(x,1),1),lagGen(x,tmin:tmax)];
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
%         tinds = tlims{i};
        tinds = usetinds(tlims{i},fs,minlen);
    else
        tinds = usetinds(tlims,fs,minlen);
    end
    x{i} = x{i}(tinds,:);
    y{i} = y{i}(tinds,:);
end
fprintf('Completed in %.3f s\n',toc(std_tm));

fprintf('Computing model prediction...\n');
% Calculate prediction
% model = [c;reshape(model,size(model,1)*size(model,2),size(model,3))];
model = reshape(model,size(model,1)*size(model,2),size(model,3));
% getting variables ready to store each trial
pred = cell(length(x));
r = NaN(length(x),size(model,3));
p = NaN(length(x),size(model,3));
rmse = NaN(length(x),size(model,3));
for i = 1:length(x)
    pred{i} = x{i}*model;
    % Calculate accuracy
    if ~isempty(y{i})
        [all_r,all_p] = corr(y{i},pred{i});
        r(i,:) = diag(all_r); p(i,:) = diag(all_p);
        rmse(i,:) = sqrt(mean((y{i}-pred{i}).^2,1))';
    end
end

end