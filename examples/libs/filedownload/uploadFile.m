function [wasSuccess, response] = uploadFile(strLocalFilename, strURLFilename, options)
%uploadFile Upload a file to web while displaying progress.
%
%   uploadFile(strLocalFilename, strURLFilename) uploads the file
%   specified by the local path `strLocalFilename` to the web location 
%   specified by `strURLFilename`.
%
%   wasSuccess = uploadFile(localFilename, strURLFilename) uploads the file
%   and returns a boolean value indicating if the upload was successful or
%   not.
%
%   [wasSuccess, response] = uploadFile(localFilename, strURLFilename) 
%   uploads the file and returns the wasSuccess boolean and a response
%   object.
%
%   Options for the progress display:
%       DisplayMode     : Where to display progress. Options: 'Dialog Box' (default) or 'Command Window'
%       UpdateInterval  : Interval (in seconds) for updating progress. Default = 1 second.
%       ShowFilename    : Whether to show name of uploaded file. Default = false.
%       IndentSize      : Size of indentation if displaying progress in command window.

%   Written by Eivind Hennestad | v1.0.6

    arguments 
        strLocalFilename       char         {mustBeNonempty}
        strURLFilename         char         {mustBeValidUrl}
        options.DisplayMode    char         {mustBeValidDisplay} = 'Dialog Box'
        options.UpdateInterval (1,1) double {mustBePositive}     = 1
        options.ShowFilename   (1,1) logical                     = false
        options.IndentSize     (1,1) uint8                       = 0
    end

    if options.ShowFilename
        [~, filename, ext] = fileparts(strURLFilename);
        filename = [char(filename), char(ext)];
    else
        filename = '';
    end

    monitorOpts = {...
        'DisplayMode', options.DisplayMode, ...
        'UpdateInterval', options.UpdateInterval, ...
        'Filename', filename, ...
        'IndentSize', options.IndentSize };
    
    webOpts = matlab.net.http.HTTPOptions(...
        'ProgressMonitorFcn', @(opts) FileTransferProgressMonitor(monitorOpts{:}),...
        'UseProgressMonitor', true, ...
        'ConnectTimeout', 20);

    % Create a file provider for uploading the file
    provider = matlab.net.http.io.FileProvider(strLocalFilename);
    
    method = matlab.net.http.RequestMethod.PUT;
    req = matlab.net.http.RequestMessage(method, [], provider);
    
    strURLFilename = matlab.net.URI(strURLFilename);
    
    [response, ~, ~] = req.send(strURLFilename, webOpts);
    
    if response.StatusCode == matlab.net.http.StatusCode.OK
        wasSuccess = true;
    else
        wasSuccess = false;
    end
    
    if nargout < 1
        if ~wasSuccess
            error(string(response.StatusLine))
        end
        clear wasSuccess
    end

    if nargout < 2
        clear response
    end
end


%% Custom validation functions

function mustBeValidDisplay(displayName)
    mustBeMember(displayName, {'Dialog Box', 'Command Window'})
end

function mustBeValidUrl(urlString)
    try
        matlab.internal.webservices.urlencode(urlString);
    catch ME
        throwAsCaller(ME)
    end
end
