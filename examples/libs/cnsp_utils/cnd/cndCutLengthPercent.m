function cnd = cndCutLengthPercent(cnd,percent2keep)
%CNDCUTSAMELENGTH cuts cnd by keeping the first portion of signal specified
%in percentage by the parameter percent2keep.
%   [CND] = CNDCUTLENGTHPERCENT(CND,PERCENT)
%
%   Author: Giovanni Di Liberto
%   Last update: 26 January 2023
%   Copyright 2023 Di Liberto Lab, Trinity College Dublin

    percent2keep = min(100,percent2keep);
    percent2keep = max(0,percent2keep);
    for tr = 1:length(cnd.data)
        len2keep = round(size(cnd.data{1,tr},1)*percent2keep/100);
        for iFea = size(cnd.data,1)
            cnd.data{iFea,tr} = cnd.data{iFea,tr}(1:len2keep,:);
        end
    end
end
    