function cnd = cndReref(cnd,type)
%CNDREREF  Re-referencing of neural data (e.g., EEG) in CND format.
%   CND = CNDREREF(CND,TYPE) re-referencing of the neural data contained in
%   the variable CND. The re-referencing can be to one particular set of
%   external channels (e.g., 'Mastoids') or to the global average ('Avg')
%   of all channels in the 'eeg.data' field. 
%
%   CNDREREF returns the CND structure after re-referencing the data
%       'cnd'       -- neural data in the Continuous-event Neural Data
%                      format (CND)
%       'type'      -- either 'Avg' or the name of the external electrodes
%                      that should be used for the referencing 'Mastoids'.
%                      If 'type' is a number, the referencing will be
%                      applied based on the corresponding channel in the 
%                      'eeg.data' field.
%
%   Author: Giovanni Di Liberto
%   Last update: 11 June 2021
%   Copyright 2021 Di Liberto Lab, Trinity College Dublin

    if isempty(cnd) || isempty(cnd.data)
        disp('The CND structure is empty or not a cell array')
        return
    elseif ~iscell(cnd.data)
        disp('The CND.data structure is not a cell array')
        return
    end
    
    % Validate input parameters
    validateparamin(cnd,type)
    
    % Checking that the re-referencing was not run already
    if isfield(cnd,'reRef')
        disp('Warning: The re-referencing was run already. Skipping this execution')
        return
    end

    % Search for re-reference channels
    if isstring(type) || ischar(type)
        if strcmp(type,'Avg') || strcmp(type,'Average')
            refVec = cell(length(cnd.data),1); % preallocate the cell array for the reference
            for tr = 1:length(cnd.data)
                refVec{tr} = mean(cnd.data{tr},2); % avg of all channels
                cnd.data{tr} = cnd.data{tr} - refVec{tr};
            end
        else
            extIdx = 0;
            extFound = 0;
            if ~isfield(cnd,'extChan')
                error('The external channel field is missing')
                return
            end

            for ii = 1:length(cnd.extChan)
                if strcmp(cnd.extChan{ii}.description,type)
                    extIdx = ii;
                    extFound = extFound + 1;
                end
            end

            if ~extIdx
                error('The external channel field was not found')
                return
            end

            if extFound > 1
                disp("Warning: the '"+type+"' description was found in correspondence of multiple external channel fields")
            end

            refVec = cnd.extChan{extIdx}.data;
            
            % Re-referencing to the reference vector (average across channels if
            % multiple channels (e.g., left and right mastoids)
            for tr = 1:length(refVec)
                cnd.data{tr} = cnd.data{tr} - mean(refVec{tr},2);
            end
            
%             % Re-referencing the external channels too
%             for ee = 1:length(cnd.extChan)
%                 for tr = 1:length(refVec)
%                     cnd.extChan{extIdx}.data{tr} = cnd.extChan{extIdx}.data{tr} - mean(refVec{tr},2);
%                 end
%             end
        end
    else % if the user input the 'type' wrong, or if none was specified
        error('Type must be either Avg or the label for external channels (like Mastoids)')
%         % Re-referencing to the reference vector (average across channels if
%         % multiple channels (e.g., Cz and FCz)
%         for tr = 1:length(cnd.data)
%             cnd.data{tr} = cnd.data{tr} - mean(cnd.data{tr}(:,type),2);
%         end
    end
    cnd.reRef = type;
    cnd = cndNewOp(cnd,"Re-referencing to "+type);
end


function validateparamin(cnd,type)
%VALIDATEPARAMIN  Validate input parameters.
%   VALIDATEPARAMIN(CND,TYPE) validates the input parameters
%   of the main function.

    nChans = cndGetNChans(cnd);
    if (isnumeric(type) || isscalar(type)) && (type <= 0 || type > nChans)
        error('When numeric, "type" must indicate a valid channel in the neural data.')
    elseif (isnumeric(type) || isscalar(type)) && nChans == -1
        error('Different number of channels for distinct subjects. Re-referencing to a particular channel is not supported')
    end
    
    if isfield(cnd,'extChan')
        for ii = 1:length(cnd.extChan)
            if length(cnd.data) ~= length(cnd.extChan{ii}.data)
                error('External channels and main data have different number of elements (trials or runs)')
            end
        end
    end
end
