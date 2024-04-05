function [vols statmat_spline] = FEMA_convert_splinevols(tbl_bf, xvec, statmat, mask, ...
                                        outpath)
% Function to convert the test statistics of a FEMA analysis using spline
% basis functions to a 4D vols variable.
%% Inputs:
% tbl_bf:       table of basis functions (typically loaded using readtable)
%
% xvec:         1-D array of x-values, length is equal to number of 
%               observations in model
%
% statmat:      array of test statistics (e.g. betas, z scores) for basis 
%               function variables, derived from FEMA output; number of 
%               columns is equal to number of basis functions used in FEMA 
%               analysis
%
% mask:         mask for use in calculating volume
%
% xvec:         vector of x values over which basis functions are evaluated;
%               length is equal to number of rows in tbl_bf
%
% outpath:      optional; path to location to save nifti file; if empty, 
%               file will not be saved
%
%% Output:
% vols:         4D array of volumes to be read by a 3D volume viewer; 4th 
%               dimension represents different values of x over which spline 
%               function is evaluated

%% calculate derivatives of basis functions
bfmat = table2array(tbl_bf);
[dummy dbfmat] = gradient(bfmat); dbfmat = dbfmat/(xvec(2)-xvec(1));

%% compute spline function as weighted sum of betas for basis functions
vols = NaN([size(mask) length(xvec)]);
beta_hat_spline = NaN([length(xvec) size(statmat,2)]);
nframes = size(statmat,2)/length(find(mask>=0.5));
if floor(nframes)~=nframes
    warning('Size of stat matrix is incompatible with running fullvol.');
end
for agei = 1:length(xvec)
    wvec = dbfmat(agei,:)';
    valvec = sum(statmat.*wvec,1);
    statmat_spline(agei, :) = valvec;
    if floor(nframes)==nframes vols(:,:,:,agei) = fullvol(valvec,mask); end
end

% if outpath exists
if exist('outpath', 'var')
    % Taken from writeNIFTI.m function
    % These should be saved along with volinfo, and passed along to this function; hardcode for now
    M_atl_sub = [0    -2     0   102; 0     0     2  -132; -2     0     0   102; 0     0     0     1]; 

    % keyboard;

    % save vols as nifti
    niftiwrite_amd(vols,outpath,M_atl_sub);
end