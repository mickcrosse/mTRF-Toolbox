function xpc = dili_ccaDataPrep_shifts(xx,NKEEP,shifts)
    for iTrial=1:length(xx)
        x = xx{iTrial};
        xx{iTrial} = nt_multishift(x,shifts);
    end
    % PCA of components and time-lags (or filter-bank)
    topcs = nt_pcarot(nt_cov(xx));
    xpc   = nt_mmat(xx,topcs(:,1:min(NKEEP,size(topcs,2)))); % reduce dimensionality
end