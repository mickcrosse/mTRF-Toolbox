function [adi,d] = mTRFattnevaluate(acc1,acc2,varargin)
%MTRFATTNEVALUATE  Evaluate the accuracy of an attention decoder.
%   ADI = MTRFATTNEVALUATE(ACC1,ACC2) returns the attention decoding
%   index (ADI) based on the proportion of folds where the accuracy for
%   attended stimuli ACC1 was greater than that of the unattended stimuli
%   ACC2.
%
%   If ACC1 or ACC2 are matrices, it is assumed that the rows correspond to
%   observations and the columns to variables, unless otherwise stated via
%   the 'dim' parameter (see below). If they are vectors, it is assumed
%   that the first non-singleton dimension corresponds to observations.
%   ACC1 and ACC2 must have the same number of observations.

%   [ADI,D] = MTRFATTNEVALUATE(ACC1,ACC2) returns the sensitivity index
%   based on d', where the accuracy for attended stimuli is considered
%   signal and that of the unattended stimuli is considered noise.
%
%   [...] = MTRFATTNEVALUATE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both ACC1 and ACC2.
%
%   See also DPRIME.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   References:
%      [1] Crosse MC, Di Liberto GM, Bednar A, Lalor EC (2016) The
%          multivariate temporal response function (mTRF) toolbox: a MATLAB
%          toolbox for relating neural signals to continuous stimuli. Front
%          Hum Neurosci 10:604.

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

% Parse input arguments
arg = parsevarargin(varargin);

% Get dimensions
if arg.dim == 2
    acc1 = acc1';
    acc2 = acc2';
end
n = size(acc1,1);
if size(acc2,1) ~= n
    error(['ACC1 and ACC2 arguments must have the same number of '...
        'observations.'])
end

% Compute means
m1 = sum(acc1,1)/n;
m2 = sum(acc2,1)/n;

% Compute ADI and d'
adi = sum(acc1>acc2,1)/n;
d = (m1 - m2)./sqrt((sum((acc1-m1).^2) + sum((acc2-m2).^2))/(n-1)/2);

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