function [Cxx,Cxy,folds] = olscovmat(x,y,lags,type,zeropad,verbose)
%OLSCOVMAT  Covariance matrices for ordinary least squares estimation.
%   [CXX,CXY] = OLSCOVMAT(X,Y,LAGS) returns the covariance matrices for
%   ordinary least squares (OLS) estimation using time-lagged features of
%   X. X and Y are matrices or cell arrays containing corresponding trials
%   of continuous data. LAGS is a vector of time lags in samples.
%
%   If X or Y are matrices, it is assumed that the rows correspond to
%   observations and the columns to variables. If they are cell arrays
%   containing multiple trials, the covariance matrices of each trial are
%   summed to produce CXX and CXY.
%
%   [CXX,CXY,FOLDS] = OLSCOVMAT(...) returns cell arrays containing the
%   individual folds in a structure with the following fields:
%       'xlag'      -- design matrices containing time-lagged features of X
%       'Cxx'       -- autocovariance matrices of XLAG
%       'Cxy'       -- cross-covariance matrices of XLAG and Y
%
%   [...] = OLSCOVMAT(X,Y,LAGS,TYPE) specifies the type of model that the
%   covariance matrices will be used to fit. Pass in 'multi' for TYPE to
%   use all lags simultaneously (default), or 'single' to use each lag
%   individually.
%
%   [...] = OLSCOVMAT(X,Y,LAGS,TYPE,ZEROPAD) specifies whether to zero-pad
%   the outer rows of the design matrix or delete them. Pass in 1 for
%   ZEROPAD to zero-pad them (default), or 0 to delete them.
%
%   [...] = OLSCOVMAT(X,Y,LAGS,TYPE,ZEROPAD,VERBOSE) specifies whether to
%   execute in verbose mode. Pass in 1 for VERBOSE to execute in verbose
%   mode (default), or 0 to execute in non-verbose mode.
%
%   See also COV, LSCOV, MLSCOVMAT.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <crossemj@tcd.ie>
%            Nate Zuk <zukn@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Set default values
if nargin < 4 || isempty(type)
    type = 'multi';
end
if nargin < 5 || isempty(zeropad)
    zeropad = true;
end
if nargin < 6 || isempty(verbose)
    verbose = true;
end

% Get dimensions
xvar = size(x{1},2);
yvar = size(y{1},2);
nfold = numel(x);
switch type
    case 'multi'
        nvar = xvar*numel(lags)+1;
        nlag = 1;
    case 'single'
        nvar = xvar+1;
        nlag = numel(lags);
end
if nargout > 2
    ncell = nfold;
else
    ncell = 1;
end

% Verbose mode
if verbose
    v = verbosemode([],0,nfold);
end

% Initialize variables
CxxInit = zeros(nvar,nvar,nlag);
CxyInit = zeros(nvar,yvar,nlag);
Cxx = CxxInit;
Cxy = CxyInit;
Cxxi = cell(ncell,1);
Cxyi = cell(ncell,1);
xlag = cell(ncell,1);
ii = 1;

for i = 1:nfold
    
    % Generate design matrix
    xlag{ii} = lagGen(x{i},lags,zeropad,1);
    
    switch type
        
        case 'multi'
            
            % Compute covariance matrices
            Cxxi{ii} = xlag{ii}'*xlag{ii};
            Cxyi{ii} = xlag{ii}'*y{i};
            
        case 'single'
            
            % Initialize cells
            Cxxi{ii} = CxxInit;
            Cxyi{ii} = CxyInit;
            
            for j = 1:nlag
                
                % Index lag
                idx = [1,xvar*(j-1)+2:xvar*j+1];
                
                % Compute covariance matrices
                Cxxi{ii}(:,:,j) = xlag{ii}(:,idx)'*xlag{ii}(:,idx);
                Cxyi{ii}(:,:,j) = xlag{ii}(:,idx)'*y{i};
                
            end
            
    end
    
    % Sum covariance matrices
    Cxx = Cxx + Cxxi{ii};
    Cxy = Cxy + Cxyi{ii};
    
    if nargout > 2
        ii = ii+1;
    end
    
    % Verbose mode
    if verbose
        v = verbosemode(v,i,nfold);
    end
    
end

% Format output
if nargout > 2
    folds = struct('xlag',{xlag},'Cxx',{Cxxi},'Cxy',{Cxyi});
end

function v = verbosemode(v,fold,nfold)
%VERBOSEMODE  Execute verbose mode.
%   V = VERBOSEMODE(V,FOLD,NFOLD) prints details about the progress of the
%   main function to the screen.

if fold == 0
    v = struct('msg',[],'h',[],'tocs',0);
    fprintf('Computing covariance matrices\n')
    v.msg = ['%d/%d [',repmat(' ',1,nfold),']\n'];
    v.h = fprintf(v.msg,fold,nfold);
else
    if fold == 1 && toc < 0.1
        pause(0.1)
    end
    v.tocs = v.tocs + toc;
    fprintf(repmat('\b',1,v.h))
    v.msg = ['%d/%d [',repmat('=',1,fold),repmat(' ',1,nfold-fold),'] - ',...
        '%.3fs/fold\n'];
    v.h = fprintf(v.msg,fold,nfold,v.tocs/fold);
end
if fold < nfold
    tic
end