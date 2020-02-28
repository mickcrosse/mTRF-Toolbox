function [Cxx,Cxy1,Cxy2,X] = mlscovmat(x,y1,y2,lags,type,split,zeropad,sumcov)
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
%   See also LSCOV, OLSCOVMAT.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <mickcrosse@gmail.com>
%            Nate Zuk <zukn@tcd.ie>
%   Copyright 2014-2020 Lalor Lab, Trinity College Dublin.

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
    X = cell(nfold,1);
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
            
            % Generate design matrix
            [X,idx] = lagGen(x{i}(iseg,:),lags,zeropad);
            X = [ones(numel(idx),1),X]; %#ok<*AGROW>
            
            % Compute covariance matrices
            switch type
                case 'multi'
                    Cxx = Cxx + X'*X;
                    Cxy1 = Cxy1 + X'*y1{i}(iseg(idx),:);
                    Cxy2 = Cxy2 + X'*y2{i}(iseg(idx),:);
                case 'single'
                    for k = 1:nlag
                        ilag = [1,xvar*(k-1)+2:xvar*k+1];
                        Cxx(:,:,k) = Cxx(:,:,k) + ...
                            X(:,ilag)'*X(:,ilag);
                        Cxy1(:,:,k) = Cxy1(:,:,k) + ...
                            X(:,ilag)'*y1{i}(iseg(idx),:);
                        Cxy2(:,:,k) = Cxy2(:,:,k) + ...
                            X(:,ilag)'*y2{i}(iseg(idx),:);
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
            
            % Generate design matrix
            [X{n},idx] = lagGen(x{i}(iseg,:),lags,zeropad);
            X{n} = [ones(numel(idx),1),X{n}];
            
            % Compute covariance matrices
            switch type
                case 'multi'
                    Cxx{n} = X{n}'*X{n}; 
                    Cxy1{n} = X{n}'*y1{i}(iseg(idx),:);
                    Cxy2{n} = X{n}'*y2{i}(iseg(idx),:);
                case 'single'
                    Cxx{n} = zeros(nvar+1,nvar+1,nlag);
                    Cxy1{n} = zeros(nvar+1,y1var,nlag);
                    Cxy2{n} = zeros(nvar+1,y2var,nlag);
                    for k = 1:nlag
                        ilag = [1,xvar*(k-1)+2:xvar*k+1];
                        Cxx{n}(:,:,k) = X{n}(:,ilag)'*X{n}(:,ilag);
                        Cxy1{n}(:,:,k) = X{n}(:,ilag)'*y1{i}(iseg(idx),:);
                        Cxy2{n}(:,:,k) = X{n}(:,ilag)'*y2{i}(iseg(idx),:);
                    end
            end
            
        end
        
    end
    
end