function [stim,resp] = cndCutSameLength(stim,resp)
%CNDCUTSAMELENGTH cuts stim and resp, two CND structures, so that their
%data fields have the same number of samples
%   [CND_STIM,CND_RESP] = CNDCUTSAMELENGTH(CND,TYPE) downsamp
%
%   CNDREREF returns the CND structure after re-referencing the data
%       'stim'       -- stim data in the Continuous-event Neural Data
%                      format (CND)
%       'resp'       -- neural data in the Continuous-event Neural Data
%                      format (CND)
%
%   Author: Giovanni Di Liberto
%   Last update: 11 July 2022
%   Copyright 2022 Di Liberto Lab, Trinity College Dublin

    for tr = 1:length(resp.data)
        minLen = min(size(stim.data{1,tr},1),size(resp.data{tr},1));
        resp.data{tr} = resp.data{tr}(1:minLen,:);
        for iFea = size(stim.data,1)
            stim.data{iFea,tr} = stim.data{iFea,tr}(1:minLen,:);
        end
    end
end
    