function [y,topcs] = nt_pcaNorm(x,NPCS)
% This function runs a PCA analysis on the matrix 'x' (neural data matrix 
% with size samples x channels). 'y' is the rotated data, 'topcs' is the
% rotation matrix.
%
% Packaged by Giovanni Di Liberto based on code from Alain de Cheveigne
% Last update: 24 June 2022
%
    % Normalize each column
    tmp=sqrt(mean(nt_unfold(x).^2)); % norm
    NN=1./tmp; 
    NN(find(tmp==0))=0;
    tmp=[];
    y=bsxfun(@times,x,NN);

    % PCA, truncate
    topcs=nt_pca0(y);
    y=nt_mmat(y,topcs(:,1:NPCS));
end