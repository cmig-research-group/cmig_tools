function [vols_dbf, varargout] = FEMA_convert_splinevols(tbl_bf, xvec, betamat, mask, coeffCovar, outpath)

% Function to convert the test statistics of a FEMA analysis using spline
% basis functions to a 4D vols variable.
%% Inputs:
% tbl_bf:       table of basis functions (typically loaded using readtable)
%
% xvec:         1-D array of x-values, length is equal to number of 
%               observations in model
%
% betamat:      array of beta estimates for basis 
%               function variables, derived from FEMA output; number of 
%               columns is equal to number of basis functions used in FEMA 
%               analysis
%
% mask:         mask for use in calculating volume
%
% xvec:         vector of x values over which basis functions are evaluated;
%               length is equal to number of rows in tbl_bf
%
% coeffCovar:   estimated coefficient covariance matrix
%
% outpath:      optional; path to location to save nifti file; if empty, 
%               file will not be saved
%
%% Outputs:
% vols_dbf:     4D array of volumes to be read by a 3D volume viewer; 4th 
%               dimension represents different values of x over which spline 
%               function is evaluated; volumes corresponding to the derivatives 
%               of spline estimates for different values of xvec
%
% sevols_dbf:   4D array of volumes corresponding to the standard error 
%               of the derivative of the spline estimates for each value of xvec
%
% zvols_dbf:    4D array of volumes corresponding to the z-scores of the 
%               derivative of the spline estimates
%
% volmat_dbf:   matrix containing the weighted sum of the betas for the 
%               derivatives of the splines
%
% vols_bf:      4D array of volumes corresponding to the spline estimates 
%               (not their derivatives) for different values of xvec
%
% sevols_bf:    4D array of volumes corresponding to the standard error of 
%               the spline function itself
%
% zvols_bf:     4D array of volumes corresponding to the z-scores of the spline 
%               estimates
%
% volmat_bf:    matrix containing the weighted sum of the betas for the 
%               spline basis functions

%% calculate derivatives of basis functions
bfmat = table2array(tbl_bf);

colnums_bf = 1:size(bfmat,2); % Identify X columns / beta coeffcients corresponding to BFs -- this needs to be done in more elegant way

statmat = betamat(colnums_bf,:); % Select betas for columns of interest -- this shouldn't be needed

[dummy dbfmat] = gradient(bfmat); dbfmat = dbfmat/(xvec(2)-xvec(1));

%% compute spline function as weighted sum of betas for basis functions
% vols_dbf = NaN([size(mask) length(xvec)]);
% sevols_dbf = NaN([size(mask) length(xvec)]);
% zvols_dbf = NaN([size(mask) length(xvec)]);
% beta_hat_spline = NaN([length(xvec) size(statmat,2)]);
nframes = size(statmat,2)/length(find(mask>=0.5));
if floor(nframes)~=nframes
    warning('Size of stat matrix is incompatible with running fullvol.');
end

% gradient of basis functions
for xi = 1:length(xvec)
    wvec = dbfmat(xi,:)';
    valvec_dbf = sum(statmat.*wvec,1); 
    sevec_dbf = rowvec(sqrt(pagemtimes(pagemtimes(wvec',coeffCovar(colnums_bf,colnums_bf,:)),wvec)));
    zvec_dbf = valvec_dbf ./ sevec_dbf;
    volmat_dbf(xi, :) = valvec_dbf;
    semat_dbf(xi,:) = sevec_dbf;
    zmat_dbf(xi, :) = zvec_dbf;

    if floor(nframes)==nframes
        vols_dbf(:,:,:,xi) = fullvol(valvec_dbf,mask); 
        sevols_dbf(:,:,:,xi) = fullvol(sevec_dbf,mask); 
        zvols_dbf(:,:,:,xi) = fullvol(zvec_dbf, mask);
    end
end

% basis functions
for xi = 1:length(xvec)
    wvec = bfmat(xi,:)';
    valvec_bf = sum(statmat.*wvec,1); 
    sevec_bf = rowvec(sqrt(pagemtimes(pagemtimes(wvec',coeffCovar(colnums_bf,colnums_bf,:)),wvec)));
    zvec_bf = valvec_bf ./ sevec_bf;
    volmat_bf(xi, :) = valvec_bf;
    semat_bf(xi,:) = sevec_bf;
    zmat_bf(xi, :) = zvec_bf;

    if floor(nframes)==nframes
        vols_bf(:,:,:,xi) = fullvol(valvec_bf,mask); 
        sevols_bf(:,:,:,xi) = fullvol(sevec_bf,mask); 
        zvols_bf(:,:,:,xi) = fullvol(zvec_bf, mask);
    end
end

if nargout>1
    varargout{1} = sevols_dbf;
    varargout{2} = zvols_dbf;
    varargout{3} = volmat_dbf;
    varargout{4} = vols_bf;
    varargout{5} = sevols_bf;
    varargout{6} = zvols_bf;
    varargout{7} = volmat_bf;
end

% figure(1); clf;
% subplot(2,2,1); imagesc(mumat_spline); colormap(blueblackred); colorbar;
% subplot(2,2,2); imagesc(semat_spline); colorbar;
% subplot(2,2,3); imagesc(mumat_spline./semat_spline); colorbar;

% if outpath exists
if exist('outpath', 'var')
    % Taken from writeNIFTI.m function
    % These should be saved along with volinfo, and passed along to this function; hardcode for now
    M_atl_sub = [0    -2     0   102; 0     0     2  -132; -2     0     0   102; 0     0     0     1]; 

    % keyboard;

    % save vols as nifti
%    niftiwrite_amd(vols,outpath,M_atl_sub); % Need to update this to write the right information to the right place
end