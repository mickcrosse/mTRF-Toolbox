function avgModel = mTRFmodelAvg(models,normFlag)
% This function returns an mTRF model structure containing the average
% TRF weights across multiple models. This function can be useful to
% average the TRF models across participants or sessions.
%
% Input:
%  models   - array of mTRF model structures (output of the mTRFtrain
%             function)
%  normFlag - if 1, the model weights are normalised before averaging
%             across models
%
% Output:
%  avgModel  - mTRF model structure with average mTRF model weights
% 
% Author: Giovanni Di Liberto
% Last update: 9 July 2021
%
    if nargin < 1 || isempty(models)
        disp("Invalid input parameters ('model' is empty)")
        return
    end
    if nargin < 2 || isempty(normFlag)
        normFlag = 0;
    end
    
    wAvg = models(1).w;
    if normFlag, wAvg = wAvg/std(wAvg(:)); end
    for sub = 2:length(models)
        w = models(sub).w;
        if normFlag, w = w/std(w(:)); end
        wAvg = wAvg + w;
    end
    wAvg = wAvg/length(models);
    
    avgModel = models(1);
    avgModel.w = wAvg;
end
