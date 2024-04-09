function [acc,d] = mTRFattnevaluate(r1,r2,varargin)
%MTRFATTNEVALUATE  Evaluate the performance of an attention decoder.
%   ACC = MTRFATTNEVALUATE(R1,R2) returns the accuracy of an attention
%   decoder based on the proportion of observations where the correlation
%   for the attended stimulus R1 was greater than that of the unattended
%   stimulus R2 as per O'Sullivan et al. (2015). Note, metrics other than
%   correlation can be input to this function to evaluate performance.
%
%   If R1 or R2 are matrices, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If they are vectors, it is assumed
%   that the first non-singleton dimension corresponds to observations.
%   R1 and R2 must have the same number of observations.
%
%   [ACC,D] = MTRFATTNEVALUATE(R1,R2) returns the attention modulation
%   index based on d', where R1 is considered signal and R2 is considered
%   noise as per de Cheveigné et al. (2018).
%
%   [...] = MTRFATTNEVALUATE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both R1 and R2.
%
%   See also DPRIME.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] O'Sullivan JA, Power AJ, Mesgarani N, Rajaram S, Foxe JJ, Shinn-
%          Cunningham BG, Slaney M, Shamma SA, Lalor EC (2015) Attentional
%          Selection in a Cocktail Party Environment Can Be Decoded from
%          Single-Trial EEG. Cereb Cortex 25(7):1697-1706.
%      [2] de Cheveigné A, Wong DE, Di Liberto GM, Hjortkjær J, Slaney M,
%          Lalor EC (2015) Decoding the auditory brain with canonical
%          component analysis. NeuroImage 172:206-216.

%   Authors: Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Get dimensions
if arg.dim == 2
    r1 = r1';
    r2 = r2';
end
n = size(r1,1);
if size(r2,1) ~= n
    error('R1 and R2 arguments must have the same number of observations.')
end

% Compute accuracy
acc = squeeze(sum(r1>r2,1)/n);

if nargout > 1
    
    % Get mean values
    m1 = sum(r1,1)/n;
    m2 = sum(r2,1)/n;
    
    % Compute d'
    d = squeeze((m1 - m2)./sqrt((sum((r1 - m1).^2) + ...
        sum((r2 - m2).^2))/(n - 1)/2));
    
end

function arg = parsevarargin(varargin)
%PARSEVARARGIN  Parse input arguments.
%   [PARAM1,PARAM2,...] = PARSEVARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   parses the input arguments of the main function.

% Create parser object
p = inputParser;

% Dimension to work along
errorMsg = 'It must be a positive integer scalar within indexing range.';
validFcn = @(x) assert(x==1||x==2,errorMsg);
addParameter(p,'dim',1,validFcn);

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;