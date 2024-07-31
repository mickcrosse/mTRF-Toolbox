function cnd = cndNewOp(cnd,operationStr)
%CNDNEWOP  Updating the CND structure by adding the name of the operation
%   that is going to be exectuted on the data. This will help keeping track
%   of the preprocessing and processing run on the data.
%
%   [] = CNDNEWOP(CND,OPERATIONSTR)
%
%   CNDNEWOP prints the operation name on the terminal
%       'cnd'          -- neural data in the Continuous-event Neural Data
%                         format (CND)
%       'operationStr' -- name of the operation to remember
%
%   Author: Giovanni Di Liberto
%   Last update: 17 June 2021
%   Copyright 2021 Di Liberto Lab, Trinity College Dublin

    if isempty(cnd)
        disp('Error: The CND structure is empty')
        return
    end

    if ~isfield(cnd,'processingPipeline')
        cnd.processingPipeline{1} = operationStr; 
    else
        cnd.processingPipeline{end+1} = operationStr; 
    end
    disp("Processing CND: "+operationStr)
end