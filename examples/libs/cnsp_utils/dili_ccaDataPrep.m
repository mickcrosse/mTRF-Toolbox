function [xpc,topcs] = dili_ccaDataPrep(xx,NKEEP,smoothBases)
    
    if ~exist('smoothBases') || isempty(smoothBases)
%         smoothBases = round(2.^unique(0:1:8));
        smoothBases = round(2.^unique(0:0.5:7));
        % smoothBases = round(2.^unique(0:0.5:7));
    end
    
    for iTrial=1:length(xx)
        x = xx{iTrial};
        xx{iTrial} = diff(nt_multismooth(x,smoothBases),[],2); % approximate log bandwidth filterbank
%         xx{iTrial} = nt_multismooth(x,smoothBases); % approximate log bandwidth filterbank
        xx{iTrial} = xx{iTrial}(:,setdiff(1:size(xx{iTrial},2),length(smoothBases):length(smoothBases):size(xx{iTrial},2))); % if multivariate, it removes diff between different variables - smoothing is within variable
    end
    % PCA of components and time-lags (or filter-bank)
    topcs = nt_pcarot(nt_cov(xx));
    xpc   = nt_mmat(xx,topcs(:,1:min(NKEEP,size(topcs,2)))); % reduce dimensionality
end