function [surf_dbf, varargout] = FEMA_convert_splinesurf(basisSubset, Xvars, beta_hat, colnums_bf, colnum_intercept, mask, coeffCovar, outpath)

% Function to convert the test statistics of a FEMA analysis using spline
% basis functions to a 4D surf variable.
%% Inputs:
% basisSubset:   table of basis functions output by createBasisFunction (typically loaded using readtable)
%
% Xvars:        vector of x values over which basis functions are evaluated output by creasteBasisFunctions.m;
%               length is equal to number of rows in basisSubset
%
% beta_hat:      array of beta estimates from FEMA output; 
%
% colnums_bf:   column numbers of beta_hat corresponding to basis functions
%
% colnum_intercept: column number of beta_hat corresponding to intercept
%
% Xvars:         vector of x values over which basis functions are evaluated;
%               length is equal to number of rows in basisSubset
%
% mask:         binary mask of vertices to be included. 
%
% coeffCovar:   estimated coefficient covariance matrix
%
% outpath:      optional; path to location to save nifti file; if empty, 
%               file will not be saved
%
%% Outputs:
% surf_dbf:     derivatives of spline estimates for different values of Xvars
%
% semat_dbf:   standard error of the derivative of the spline estimates for each value of Xvars
%
% zmat_dbf:    z-scores of the derivative of the spline estimates
%
% surf_bf:      the spline estimates (not their derivatives) for different values of Xvars
%
% semat_bf:    the standard error of the spline function itself
%
% zmat_bf:     z-scores of the spline estimates
%


%% calculate derivatives of basis functions
bfmat = table2array(basisSubset);

if ~exist('colnums_bf', 'var')
	colnums_bf = 1:size(bfmat,2); % Identify X columns / beta coeffcients corresponding to BFs -- this needs to be done in more elegant way
end 

statmat_dbf = beta_hat(colnums_bf,:); % Select betas for columns of interest -- this shouldn't be needed
statmat_bf = beta_hat([colnums_bf colnum_intercept],:);
ivec_mask = find(mask);
[dummy dbfmat] = gradient(bfmat); 
dbfmat = dbfmat/(Xvars(2)-Xvars(1));

%% compute spline function as weighted sum of betas for basis functions

% gradient of basis functions
semat_dbf = NaN([length(Xvars) size(beta_hat,2)]);
zmat_dbf = NaN([length(Xvars) size(beta_hat,2)]);
for xi = 1:length(Xvars)
    wvec = dbfmat(xi,:)';
    valvec_dbf = sum(statmat_dbf.*wvec,1); 
    sevec_dbf = rowvec(sqrt(pagemtimes(pagemtimes(wvec',coeffCovar(colnums_bf,colnums_bf,:)),wvec)));
    zvec_dbf = valvec_dbf(ivec_mask) ./ sevec_dbf;
    surf_dbf(xi, :) = valvec_dbf;
    semat_dbf(xi, ivec_mask) = sevec_dbf;
    zmat_dbf(xi, ivec_mask) = zvec_dbf;
end

% basis functions
semat_bf = NaN([length(Xvars) size(beta_hat,2)]);
zmat_bf = NaN([length(Xvars) size(beta_hat,2)]);
for xi = 1:length(Xvars)
    wvec = bfmat(xi,:)';
    valvec_bf = sum(statmat_bf.*[wvec; 1],1); 
    sevec_bf = rowvec(sqrt(pagemtimes(pagemtimes([wvec; 1]',coeffCovar([colnums_bf colnum_intercept], [colnums_bf colnum_intercept],:)),[wvec; 1])));
    zvec_bf = valvec_bf(ivec_mask) ./ sevec_bf;
    surf_bf(xi, :) = valvec_bf;
    semat_bf(xi, ivec_mask) = sevec_bf;
    zmat_bf(xi, ivec_mask) = zvec_bf;

end

if nargout>1
    varargout{1} = semat_dbf;
    varargout{2} = zmat_dbf;
    varargout{4} = surf_bf;
    varargout{5} = semat_bf;
    varargout{6} = zmat_bf;
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

    % save surf as nifti
%    niftiwrite_amd(surf,outpath,M_atl_sub); % Need to update this to write the right information to the right place
end