function [AA,BB,RR,iBest,accMM]=nt_cca_crossvalidate_mm(xx,yy,shifts,windowSize,ncomp,A0,B0)
%[AA,BB,RR,iBest]=nt_cca_crossvalidate(xx,yy,shifts,ncomp,A0,B0) - CCA with cross-validation
%
%  AA, BB: cell arrays of transform matrices
%  RR: r scores (3D) for all components, shifts and trials
%  iBest: index of best shift
%
%  xx,yy: cell arrays of column matrices
%  shifts: array of shifts to apply to y relative to x (can be negative)
%  ncomp: number of components to consider for iBest [default: all]
%  A0,B0: if present, use these CCA transform matrices 
%
%  Plot correlation re shifts for matching trials
%    plot(shifts, mean(RR,3)');
%  Plot mean correlation re shifts for mismatched trials
%    plot(shifts, mean(mean(RR,4),3)');
%
% This version also evaluates the model according to a match-vs-mismatch
% classification.
% Modified by Giovanni M. Di Liberto
% 24 June 2022

if nargin<6
    A0=[]; B0=[]; 
end
if nargin<5; ncomp=[]; end
if nargin<3 || isempty (shifts); shifts=[0]; end
if nargin<2; error('!'); end
if ~iscell(xx) || ~iscell(yy); error('!'); end
if length(xx) ~= length (yy); error('!'); end
%if size(xx{1},1) ~= size(yy{1},1); error('!'); end
if size(xx{1},1) ~= size(yy{1},1); 
    for iTrial=1:numel(xx);
        tmp=min(size(xx{iTrial},1),size(yy{iTrial},1));
        xx{iTrial}=xx{iTrial}(1:tmp,:);
        yy{iTrial}=yy{iTrial}(1:tmp,:);
    end
end

nTrials=length(xx);

if isempty(A0)
    % calculate covariance matrices
    n=size(xx{1},2)+size(yy{1},2);
    C=zeros(n,n,length(shifts),nTrials);
    disp('Calculate all covariances...'); tic;
    nt_whoss;
    for iTrial=1:nTrials
        C(:,:,:,iTrial)=nt_cov_lags(xx{iTrial}, yy{iTrial},shifts);
    end

    % calculate leave-one-out CCAs
    disp('Calculate CCAs...'); tic;
    for iTrial=1:nTrials
        CC=sum(C(:,:,:,setdiff(1:nTrials,iTrial)),4); % covariance of all trials except iOut
        [A,B,R]=nt_cca([],[],[],CC,size(xx{1},2));  % CCA to apply to that trial (trained on others)
        AA{iTrial}=A;
        BB{iTrial}=B;
    end
    clear C CC
    toc;
else
    % set to given values
    for iTrial=1:nTrials
        AA{iTrial}=A0;
        BB{iTrial}=B0;
    end
end

%%
% calculate leave-one-out correlation coefficients
disp('Calculate cross-correlations...'); tic;
clear accMM
for iShift=1:length(shifts)
    countCorrectClassif = 0;
    countAllClassif = 0;
    
    xxx={}; yyy={};
    % shift, trim to same length, convert to CCs, normalize
    for iTrial=1:nTrials
        [xxx{iTrial},yyy{iTrial}]=nt_relshift(xx{iTrial},yy{iTrial},shifts(iShift));
        xxx{iTrial}=nt_normcol( nt_demean( nt_mmat(xxx{iTrial},AA{iTrial}(:,:,iShift)) ) );
        yyy{iTrial}=nt_normcol( nt_demean( nt_mmat(yyy{iTrial},BB{iTrial}(:,:,iShift)) ) );
    end
    for iTrial=1:nTrials
        x=xxx{iTrial}; % x in MCC space
        y=yyy{iTrial}; % y in MCC space
        % Circshift
        shiftSamples = size(x,1)/2; % circular shifts of half of the length of the trial +- 5% samples
        shiftSamples = round((shiftSamples + 2*(rand-0.5)*0.05*shiftSamples));
        xMM = x(circshift(1:size(x,1),shiftSamples),:);
        
        % Correlation
        RR(:,iShift,iTrial)=diag(x'*y) / size(x,1);
        % Match vs. mismatch classification
        % Chunking
        chunkSt = 1:windowSize:size(x,1)-windowSize;
        clear sstim xChunks yChunks xChunksMM
        for iChunk = 1:length(chunkSt)
            st = chunkSt(iChunk);
            fin = st+windowSize-1;
            xChunks(:,:,iChunk) = x(st:fin,:);
            yChunks(:,:,iChunk) = y(st:fin,:);
            xChunksMM(:,:,iChunk) = xMM(st:fin,:);
        end
%         yChunksMM = yChunks(:,:,randperm(size(yChunks,3)));
        if ~exist('ncomp') || ncomp <=0
            ncomp = size(xChunks,2);
        end
        % Classification
        for iChunk = 1:length(chunkSt)
            x1 = xChunks(:,1:ncomp,iChunk); %x1 = x1(:); % using all 1:ncomp components
            y1 = yChunks(:,1:ncomp,iChunk); %y1 = y1(:);
            x2 = xChunksMM(:,1:ncomp,iChunk); %x2 = x2(:);
            
            rMall = zeros(ncomp,1); rMMall = rMall;
            for cc = 1:ncomp
                rMall(cc) = corr(x1(:,cc),y1(:,cc));
                rMMall(cc) = corr(x2(:,cc),y1(:,cc));
            end
            if mean(rMall) > mean(rMMall) %corr(x1,y1) > corr(x2,y1)
                countCorrectClassif = countCorrectClassif + 1;
            end
        end
        countAllClassif = countAllClassif + length(chunkSt);
    end
    accMM(iShift) = countCorrectClassif/countAllClassif;
end
toc;

if isempty(ncomp); ncomp=size(RR,1); end
[~,iBest]=max(mean(mean(RR(1:ncomp,:,:),3),1)'); 

disp('done');

