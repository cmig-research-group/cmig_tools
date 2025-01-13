function [vols_dbf, varargout] = FEMA_convert_splinevols(basisSubset, Xvars, beta_hat, colnums_bf, colnum_intercept, mask, coeffCovar, outpath)

% Function to convert the test statistics of a FEMA analysis using spline
% basis functions to a 4D vols variable.
%% Inputs:
% basisSubset:       table of basis functions output by createBasisFunction (typically loaded using readtable)
%
% beta_hat:      array of beta estimates from FEMA output; 
%
% colnums_bf:   column numbers of beta_hat corresponding to basis functions
%
% colnum_intercept: column number of beta_hat corresponding to intercept
%
% mask:         mask for use in calculating volume
%
% Xvars:         vector of x values over which basis functions are evaluated output by creasteBasisFunctions.m;
%               length is equal to number of rows in basisSubset
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
%               of spline estimates for different values of Xvars
%
% sevols_dbf:   4D array of volumes corresponding to the standard error 
%               of the derivative of the spline estimates for each value of Xvars
%
% zvols_dbf:    4D array of volumes corresponding to the z-scores of the 
%               derivative of the spline estimates
%
% volmat_dbf:   matrix containing the weighted sum of the betas for the 
%               derivatives of the splines
%
% vols_bf:      4D array of volumes corresponding to the spline estimates 
%               (not their derivatives) for different values of Xvars
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
bfmat = table2array(basisSubset);

if ~exist('colnums_bf', 'var')
	colnums_bf = 1:size(bfmat,2); % Identify X columns / beta coeffcients corresponding to BFs -- this needs to be done in more elegant way
end 

statmat_dbf = beta_hat(colnums_bf,:); % Select betas for columns of interest -- this shouldn't be needed
statmat_bf = beta_hat([colnums_bf colnum_intercept],:);

[dummy dbfmat] = gradient(bfmat); 
dbfmat = dbfmat/(Xvars(2)-Xvars(1));

%% compute spline function as weighted sum of betas for basis functions
nframes = size(statmat_dbf,2)/length(find(mask>=0.5));
if floor(nframes)~=nframes
    warning('Size of stat matrix is incompatible with running fullvol.');
end

% gradient of basis functions
volmat_dbf = NaN(length(Xvars), size(statmat_dbf,2));
semat_dbf = NaN(length(Xvars), size(statmat_dbf,2));
zmat_dbf = NaN(length(Xvars), size(statmat_dbf,2));
for xi = 1:length(Xvars)
    wvec = dbfmat(xi,:)';
    valvec_dbf = sum(statmat_dbf.*wvec,1); 
    sevec_dbf = rowvec(sqrt(pagemtimes(pagemtimes(wvec',coeffCovar(colnums_bf,colnums_bf,:)),wvec)));
    zvec_dbf = valvec_dbf ./ sevec_dbf;
    volmat_dbf(xi, :) = valvec_dbf;
    semat_dbf(xi,:) = sevec_dbf;
    zmat_dbf(xi, :) = zvec_dbf;
end
vols_dbf = fullvol(volmat_dbf,mask); 
sevols_dbf = fullvol(semat_dbf,mask); 
zvols_dbf = fullvol(zmat_dbf, mask);

% basis functions
volmat_bf = NaN(length(Xvars), size(statmat_bf,2));
semat_bf = NaN(length(Xvars), size(statmat_bf,2));
zmat_bf = NaN(length(Xvars), size(statmat_bf,2));
for xi = 1:length(Xvars)
    wvec = bfmat(xi,:)';
    valvec_bf = sum(statmat_bf.*[wvec; 1],1); 
    sevec_bf = rowvec(sqrt(pagemtimes(pagemtimes([wvec; 1]',coeffCovar([colnums_bf colnum_intercept], [colnums_bf colnum_intercept],:)),[wvec; 1])));
    zvec_bf = valvec_bf ./ sevec_bf;
    volmat_bf(xi, :) = valvec_bf;
    semat_bf(xi,:) = sevec_bf;
    zmat_bf(xi, :) = zvec_bf;
end
vols_bf = fullvol(volmat_bf,mask); 
sevols_bf = fullvol(semat_bf,mask); 
zvols_bf = fullvol(zmat_bf, mask);

if nargout>1
    varargout{1} = sevols_dbf;
    varargout{2} = zvols_dbf;
    varargout{3} = volmat_dbf;
    varargout{4} = vols_bf;
    varargout{5} = sevols_bf;
    varargout{6} = zvols_bf;
    varargout{7} = volmat_bf; 
	varargout{8} = semat_bf;
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