function y = rampstim(x,Fs,varargin)
% Apply a raised sine ramp to the start and end of x

rtime = 15; % ramp time, in ms;

% Parse varargin
if ~isempty(varargin),
    for n=2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

% Make sure time is dimension 1
if size(x,1)<size(x,2), x = x'; end

% Make starting and ending ramps
tramp = (0:floor(rtime/1000*Fs)-1)/floor(rtime/1000*Fs); % time, relative to 1 cycle of cosine function
ramp = 0.5*(1-cos(pi*tramp)); % ramp shape
wnd = [ramp'; ones(length(x)-2*length(ramp),1); fliplr(ramp)'];

% Apply ramp to all channels
y = x.*(wnd*ones(1,size(x,2)));