function [Cxx,Cxy,xlag] = olscovmat(x,y,lags,type,split,zeropad,sumcov)
%OLSCOVMAT  Ordinary least squares covariance matrix estimation.
%   [CXX,CXY] = OLSCOVMAT(X,Y,LAGS) returns the covariance matrices for
%   ordinary least squares (OLS) estimation using time-lagged features of
%   X. X and Y are matrices or cell arrays containing corresponding trials
%   of continuous data, with columns corresponding to observations and rows
%   corresponding to variables. LAGS is a scalar or vector of time lags in
%   samples.
%
%   [CXX,CXY,XLAG] = OLSCOVMAT(...) returns the time-lagged features of X
%   used to compute the covariance matrices.
%
%   [...] = OLSCOVMAT(X,Y,LAGS,TYPE) specifies the type of model that the
%   covariance matrices will be used to fit. Pass in 'multi' for TYPE to
%   use all lags simultaneously (default), or 'single' to use each lag
%   individually.
%
%   [...] = OLSCOVMAT(X,Y,LAGS,TYPE,SPLIT) specifies the number of segments
%   in which to split each trial of data when computing the covariance
%   matrices. This is useful for reducing memory usage on large datasets.
%   By default, the entire trial is used.
%
%   [...] = OLSCOVMAT(X,Y,LAGS,TYPE,SPLIT,ZEROPAD) specifies whether to
%   zero-pad the outer rows of the design matrix or delete them. Pass in
%   1 to zero-pad them (default), or 0 to delete them.
%
%   [...] = OLSCOVMAT(X,Y,LAGS,TYPE,SPLIT,ZEROPAD,SUMCOV) specifies whether
%   to sum over the covariance matrices of each trial or return each one in
%   a separate cell. Pass in 1 to sum them (default), or 0 to return them
%   separately.
%
%   See also MLSCOVMAT.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse, Nate Zuk
%   Contact: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Jan 2020; Last revision: 11-Feb-2020

% Set default values
if nargin < 4 || isempty(type)
    type = 'multi';
end
if nargin < 5 || isempty(split)
    split = 1;
end
if nargin < 6 || isempty(zeropad)
    zeropad = true;
end
if nargin < 7 || isempty(sumcov)
    sumcov = true;
end

% Get dimensions
[~,xobs,xvar] = formatcells(x,1,0);
[~,~,yvar] = formatcells(y,1,0);
nlag = numel(lags);
xvar = unique(xvar);
yvar = unique(yvar);
nvar = xvar*nlag;
ntrial = numel(x);
nfold = ntrial*split;

% Initialize covariance matrices
if sumcov
    switch type
        case 'multi'
            Cxx = zeros(nvar+1,nvar+1);
            Cxy = zeros(nvar+1,yvar);
        case 'single'
            Cxx = zeros(xvar+1,xvar+1,nlag);
            Cxy = zeros(xvar+1,yvar,nlag);
    end
else
    xlag = cell(nfold,1);
    Cxx = cell(nfold,1);
    Cxy = cell(nfold,1);
end

if sumcov % sum over trials
    
    for i = 1:ntrial
        
        % Max segment size
        seg = ceil(xobs(i)/split);
        
        for j = 1:split
            
            % Segment indices
            iseg = seg*(j-1)+1:min(seg*j,xobs(i));
            
            % Generate time-lagged features
            [xlag,idx] = lagGen(x{i}(iseg,:),lags,zeropad);
            xlag = [ones(numel(idx),1),xlag];
            
            % Compute covariance matrices
            switch type
                case 'multi'
                    Cxx = Cxx + xlag'*xlag;
                    Cxy = Cxy + xlag'*y{i}(iseg(idx),:);
                case 'single'
                    for k = 1:nlag
                        ilag = [1,xvar*(k-1)+2:xvar*k+1];
                        Cxx(:,:,k) = Cxx(:,:,k) + ...
                            xlag(:,ilag)'*xlag(:,ilag);
                        Cxy(:,:,k) = Cxy(:,:,k) + ...
                            xlag(:,ilag)'*y{i}(iseg(idx),:);
                    end
            end
            
        end
        
    end
    
else % keep trials separate
    
    n = 0;
    for i = 1:ntrial
        
        % Max segment size
        seg = ceil(xobs(i)/split);
        
        for j = 1:split
            
            n = n+1;
            
            % Segment indices
            iseg = seg*(j-1)+1:min(seg*j,xobs(i));
            
            % Generate time-lagged features
            [xlag{n},idx] = lagGen(x{i}(iseg,:),lags,zeropad);
            xlag{n} = [ones(numel(idx),1),xlag{n}];
            
            % Compute covariance matrices
            switch type
                case 'multi'
                    Cxx{n} = xlag{n}'*xlag{n}; %#ok<*AGROW>
                    Cxy{n} = xlag{n}'*y{i}(iseg(idx),:);
                case 'single'
                    Cxx{n} = zeros(xvar+1,xvar+1,nlag);
                    Cxy{n} = zeros(xvar+1,yvar,nlag);
                    for k = 1:nlag
                        ilag = [1,xvar*(k-1)+2:xvar*k+1];
                        Cxx{n}(:,:,k) = xlag{n}(:,ilag)'*xlag{n}(:,ilag);
                        Cxy{n}(:,:,k) = xlag{n}(:,ilag)'*y{i}(iseg(idx),:);
                    end
            end
            
        end
        
    end
    
end