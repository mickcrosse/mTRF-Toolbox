function [Cxx,Cxy1,Cxy2,xlag] = mlscovmat(x,y1,y2,lags,type,split,zeropad,sumcov)
%MLSCOVMAT  Multisensory least squares covariance matrix estimation.
%   [CXX,CXY1,CXY2] = MLSCOVMAT(X,Y1,Y2,LAGS) returns the covariance
%   matrices for multisensory least squares (MLS) estimation using time-
%   lagged features of X. X, Y1 and Y2 are matrices or cell arrays
%   containing corresponding trials of continuous data, with columns
%   corresponding to observations and rows corresponding to variables. LAGS
%   is a scalar or vector of time lags in samples.
%
%   [CXX,CXY1,CXY2,XLAG] = MLSCOVMAT(...) returns the time-lagged features
%   of X used to compute the covariance matrices.
%
%   [...] = MLSCOVMAT(X,Y,LAGS,TYPE) specifies the type of model that the
%   covariance matrices will be used to fit. Pass in 'multi' for TYPE to
%   use all lags simultaneously (default), or 'single' to use each lag
%   individually.
%
%   [...] = MLSCOVMAT(X,Y,LAGS,TYPE,SPLIT) specifies the number of segments
%   in which to split each trial of data when computing the covariance
%   matrices. This is useful for reducing memory usage on large datasets.
%   By default, the entire trial is used.
%
%   [...] = MLSCOVMAT(X,Y,LAGS,TYPE,SPLIT,ZEROPAD) specifies whether to
%   zero-pad the outer rows of the design matrix or delete them. Pass in 1
%   to zero-pad them (default), or 0 to delete them.
%
%   [...] = MLSCOVMAT(X,Y,LAGS,TYPE,SPLIT,ZEROPAD,SUMCOV) specifies whether
%   to sum over the covariance matrices of each trial or return each one in
%   a separate cell. Pass in 1 to sum them (default), or 0 to return them
%   separately.
%
%   See also OLSCOVMAT.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse, Nate Zuk
%   Contact: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Jan 2020; Last revision: 11-Feb-2020

% Set default values
if nargin < 5 || isempty(type)
    type = 'multi';
end
if nargin < 6 || isempty(split)
    split = 1;
end
if nargin < 7 || isempty(zeropad)
    zeropad = true;
end
if nargin < 8 || isempty(sumcov)
    sumcov = true;
end

% Get dimensions
[~,xobs,xvar] = formatcells(x,1,0);
[~,~,y1var] = formatcells(y1,1,0);
[~,~,y2var] = formatcells(y1,1,0);
nlag = numel(lags);
xvar = unique(xvar);
y1var = unique(y1var);
y2var = unique(y2var);
nvar = xvar*nlag;
ntrial = numel(x);
nfold = ntrial*split;

% Initialize covariance matrices
if sumcov
    switch type
        case 'multi'
            Cxx = zeros(nvar+1,nvar+1);
            Cxy1 = zeros(nvar+1,y1var);
            Cxy2 = zeros(nvar+1,y2var);
        case 'single'
            Cxx = zeros(xvar+1,xvar+1,nlag);
            Cxy1 = zeros(xvar+1,y1var,nlag);
            Cxy2 = zeros(xvar+1,y2var,nlag);
    end
else
    xlag = cell(nfold,1);
    Cxx = cell(nfold,1);
    Cxy1 = cell(nfold,1);
    Cxy2 = cell(nfold,1);
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
                    Cxy1 = Cxy1 + xlag'*y1{i}(iseg(idx),:);
                    Cxy2 = Cxy2 + xlag'*y2{i}(iseg(idx),:);
                case 'single'
                    for k = 1:nlag
                        ilag = [1,xvar*(k-1)+2:xvar*k+1];
                        Cxx(:,:,k) = Cxx(:,:,k) + ...
                            xlag(:,ilag)'*xlag(:,ilag);
                        Cxy1(:,:,k) = Cxy1(:,:,k) + ...
                            xlag(:,ilag)'*y1{i}(iseg(idx),:);
                        Cxy2(:,:,k) = Cxy2(:,:,k) + ...
                            xlag(:,ilag)'*y2{i}(iseg(idx),:);
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
                    Cxy1{n} = xlag{n}'*y1{i}(iseg(idx),:);
                    Cxy2{n} = xlag{n}'*y2{i}(iseg(idx),:);
                case 'single'
                    Cxx{n} = zeros(nvar+1,nvar+1,nlag);
                    Cxy1{n} = zeros(nvar+1,y1var,nlag);
                    Cxy2{n} = zeros(nvar+1,y2var,nlag);
                    for k = 1:nlag
                        ilag = [1,xvar*(k-1)+2:xvar*k+1];
                        Cxx{n}(:,:,k) = xlag{n}(:,ilag)'*xlag{n}(:,ilag);
                        Cxy1{n}(:,:,k) = xlag{n}(:,ilag)'*y1{i}(iseg(idx),:);
                        Cxy2{n}(:,:,k) = xlag{n}(:,ilag)'*y2{i}(iseg(idx),:);
                    end
            end
            
        end
        
    end
    
end