function strLocalFilename = downloadFile(strLocalFilename, strURLFilename, options)
%downloadFile Download and save a file from web while displaying progress.
%
%   downloadFile(strLocalFilename, strURLFilename) downloads the file
%   specified by the url strURLFilename to the local path specified by
%   strLocalFile
%
%   strLocalFilename = downloadFile(localFilename, strURLFilename) 
%   downloads the file and returns the absolute path of the downloaded file
%
%   Options for the progress display:
%       DisplayMode     : Where to display progress. Options: 'Dialog Box' (default) or 'Command Window'
%       UpdateInterval  : Interval (in seconds) for updating progress. Default = 1 second.
%       ShowFilename    : Whether to show name of downloaded file. Default = false.
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

    % Create a file consumer for saving the file
    consumer = matlab.net.http.io.FileConsumer(strLocalFilename);
    
    method = matlab.net.http.RequestMethod.GET;
    req = matlab.net.http.RequestMessage(method, [], []);
    
    strURLFilename = matlab.net.URI(strURLFilename);
    
    [resp, ~, ~] = req.send(strURLFilename, webOpts, consumer);

    strLocalFilename = resp.Body.Data;

    if nargout < 1
        clear strLocalFilename
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
