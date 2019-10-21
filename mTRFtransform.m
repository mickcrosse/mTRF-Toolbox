function [f_model] = mTRFtransform(b_model,resp,recon,tlims,fs)
%mTRFtransform mTRF Toolbox transformation function.
%   F_MODEL = MTRFTRANSFORM(RECON,RESP,B_MODEL) tansforms the coefficients 
%   of the backward model B_MODEL to the coefficents of the corresponding
%   forward model F_MODEL as described in Haufe et al., (2014). This simple
%   procedure enables the neurophysiological interpretation of the 
%   weights of linear backward models.
%
%   Inputs:
%   b_model - linear backward model (chans by lags by feats)
%   resp    - neural response data (time by channels)
%   recon   - reconstructed stimulus estimate (time by features)
%      ** It is assumed that recon is already limited to the range of tlims
%      provided
%
%   Outputs:
%   f_model - transformed model weights (lags by chans)
%
%   See README for examples of use.
%
%   See also LAGGEN MTRFTRAIN MTRFPREDICT MTRFCROSSVAL MTRFMULTICROSSVAL.

%   References:
%      [1] Haufe S, Meinecke F, Gorgen K, Dahne S, Haynes JD, Blankertz B,
%          Bieﬂmann F (2014) On the interpretation of weight vectors of
%          linear models in multivariate neuroimaging. NeuroImage 87:96-110.

%   Authors: Adam Bednar, Emily Teoh, Giovanni Di Liberto, Mick Crosse
%   Email: mickcrosse@gmail.com, edmundlalor@gmail.com
%   Website: www.lalorlab.net
%   Lalor Lab, Trinity College Dublin, IRELAND
%   April 2014; Last revision: 22-Mar-2019

if nargin<4, 
    tlims = []; 
end

% If the response variable isn't a cell array (only one trial), turn it
% into a one element cell array
if ~iscell(resp), resp = {resp}; end
if ~iscell(recon), recon = {recon}; end

% Only retain sections of the response that are included in tlims
if ~isempty(tlims)
    for n = 1:length(resp)
        if iscell(tlims), % if tlims is a cell array, it means that specific indexes were supplied
            tinds = usetinds(tlims{i},fs,length(resp{n}));
        else
            tinds = usetinds(tlims,fs,length(resp{n}));
        end
        resp{n} = resp{n}(tinds,:);
    end
end

% Get the autocov matrix for resp and recon
res_t_res = compute_linreg_matrices(resp);
rec_t_rec = compute_linreg_matrices(recon);

% Transform model weights
% f_model = (resp'*resp)*b_model/(recon'*recon);
f_model = res_t_res*b_model/rec_t_rec;

% Format model weights
f_model = fliplr(f_model);

end