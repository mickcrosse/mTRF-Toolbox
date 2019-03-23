function [f_model] = mTRFtransform(b_model,resp,recon)
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

% Transform model weights
f_model = (resp'*resp)*b_model/(recon'*recon);

% Format model weights
f_model = fliplr(f_model);

end