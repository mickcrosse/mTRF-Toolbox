function modelsummary(model)
%MODELSUMMARY  Print model summary.
%   MODELSUMMARY(MODEL) prints a summary table of the model's parameters
%   and their shape.
%
%   mTRF-Toolbox https://github.com/mickcrosse/mTRF-Toolbox

%   Authors: Mick Crosse <crossemj@tcd.ie>
%   Copyright 2014-2024 Lalor Lab, Trinity College Dublin.

% Get dimensions
if isfield(model,'b')
    b1 = size(model.b,1); b2 = size(model.b,2); b3 = size(model.b,3);
    lb = numel(num2str(b2)) + numel(num2str(b2));
    if b3 > 1, lb = lb + numel(num2str(b3)); end
    nb = b1*b2*b3;
else
    nb = 0;
end
w1 = size(model.w,1); w2 = size(model.w,2); w3 = size(model.w,3);
lw = numel(num2str(w1)) + numel(num2str(w2));
if w3 > 1, lw = lw + numel(num2str(w3)); end
nw = w1*w2*w3;

% Print model
if model.Dir == 1
    fprintf('\nModel: "encoding"\n')
elseif model.Dir == -1
    fprintf('\nModel: "decoding"\n')
end

% Print table
for i = 1:52, fprintf('_'), end, fprintf('\n')
fprintf('Param (type)\t\tOutput Shape\t\t\tParam #\n')
for i = 1:52, fprintf('='), end, fprintf('\n')
if isfield(model,'b')
    if b3 > 1
        if lb < 6
            fprintf('bias\t\t\t\t(%d, %d, %d)\t\t\t\t%d\n',b1,b2,b3,nb)
        else
            fprintf('bias\t\t\t\t(%d, %d, %d)\t\t\t%d\n',b1,b2,b3,nb)
        end
    else
        if lb < 4
            fprintf('bias\t\t\t\t(%d, %d)\t\t\t\t\t%d\n',b1,b2,nb)
        else
            fprintf('bias\t\t\t\t(%d, %d)\t\t\t\t%d\n',b1,b2,nb)
        end
    end
    for i = 1:52, fprintf('_'), end, fprintf('\n\n')
end
if w3 > 1
    if lw < 6
        fprintf('weights\t\t\t\t(%d, %d, %d)\t\t\t\t%d\n',w1,w2,w3,nw)
    else
        fprintf('weights\t\t\t\t(%d, %d, %d)\t\t\t%d\n',w1,w2,w3,nw)
    end
else
    if lw < 4
        fprintf('weights\t\t\t\t(%d, %d)\t\t\t\t\t%d\n',w1,w2,nw)
    else
        fprintf('weights\t\t\t\t(%d, %d)\t\t\t\t%d\n',w1,w2,nw)
    end
end
for i = 1:52, fprintf('='), end, fprintf('\n')
fprintf('Trainable params: %d\n',nw+nb)
fprintf('Lag type: %s-lag\n',model.type)
for i = 1:52, fprintf('_'), end, fprintf('\n')