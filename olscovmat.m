function [Cxx,Cxy,xlag] = olscovmat(x,y,lags,type,split,zeropad,sumcov)
%OLSCOVMAT  Ordinary least squares covariance matrix estimation.
%   [CXX,CXY] = OLSCOVMAT(X,Y,LAGS) returns a the covariance matrices for
%   ordinary least squares (OLS) regression using time-lagged features of
%   X. X and Y are matrices or cell arrays containing corresponding trials 
%   of continuous training data, with columns corresponding to observations 
%   and rows corresponding to variables. LAGS is a scalar or vector of time 
%   lags in samples.
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
%   1 to zero-pad (default), or 0 to delete.
%
%   [...] = OLSCOVMAT(X,Y,LAGS,TYPE,SPLIT,ZEROPAD,SUMCOV) specifies whether
%   to sum over the covariance matrices of each trial or return each one in
%   a separate cell. Pass in 1 to sum (default), or 0 to return separately.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse, Nate Zuk
%   Contact: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Lalor Lab, Trinity College Dublin, IRELAND
%   Jan 2020; Last revision: 04-Feb-2020

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
nvar = xvar*nlag+1;
ntrial = numel(x);

% Initialize covariance matrices
if sumcov
    switch type
        case 'multi'
            Cxx = 0; Cxy = 0;
        case 'single'
            Cxx = zeros(xvar+1,xvar+1,nlag);
            Cxy = zeros(xvar+1,yvar,nlag);
    end
else
    xlag = cell(ntrial*split,1);
    Cxx = cell(ntrial*split,1);
    Cxy = cell(ntrial*split,1);
end

if sumcov % sum over trials
    
    for i = 1:ntrial
        
        % Get segment size
        nseg = ceil(xobs(i)/split);
        
        for j = 1:split % split trial into segments
            
            % Get segment indices
            iseg = nseg*(j-1)+1:min(nseg*j,xobs(i));
            
            % Generate time-lagged features
            [xlag,ilag] = lagGen(x{i}(iseg,:),lags,zeropad);
            xlag = [ones(numel(ilag),1),xlag];
            iseg = iseg(ilag);
            
            % Compute covariance matrices
            switch type
                case 'multi'
                    Cxx = Cxx + xlag'*xlag;
                    Cxy = Cxy + xlag'*y{i}(iseg,:);
                case 'single'
                    jj = 0;
                    for ii = 2:xvar:nvar
                        jj = jj+1;
                        idx = [1,ii:ii+xvar-1];
                        Cxx(:,:,jj) = Cxx(:,:,jj) + ...
                            xlag(:,idx)'*xlag(:,idx);
                        Cxy(:,:,jj) = Cxy(:,:,jj) + ...
                            xlag(:,idx)'*y{i}(iseg,:);
                    end
            end
            
        end
        
    end
    
else % keep trials separate
    
    m = 0;
    for i = 1:ntrial
        
        % Get segment size
        nseg = ceil(xobs(i)/split);
        
        for j = 1:split % split trial into segments
            
            m = m+1;
            
            % Get segment indices
            iseg = nseg*(j-1)+1:min(nseg*j,xobs(i));
            
            % Generate time-lagged features
            [xlag{m},ilag] = lagGen(x{i}(iseg,:),lags,zeropad);
            xlag{m} = [ones(numel(ilag),1),xlag{m}]; 
            iseg = iseg(ilag);
            
            % Compute covariance matrices
            switch type
                case 'multi'
                    Cxx{m} = xlag{m}'*xlag{m}; %#ok<*AGROW>
                    Cxy{m} = xlag{m}'*y{i}(iseg,:);
                case 'single'
                    jj = 0;
                    Cxx{m} = zeros(xvar+1,xvar+1,nlag);
                    Cxy{m} = zeros(xvar+1,yvar,nlag);
                    for ii = 2:xvar:nvar
                        jj = jj+1;
                        idx = [1,ii:ii+xvar-1];
                        Cxx{m}(:,:,jj) = xlag{m}(:,idx)'*xlag{m}(:,idx);
                        Cxy{m}(:,:,jj) = xlag{m}(:,idx)'*y{i}(iseg,:);
                    end
            end
            
        end
        
    end
    
end