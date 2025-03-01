function [beta_bf, SE_bf, Z_bf, beta_dbf, SE_dbf, Z_dbf] = FEMA_convert_splines(basisSubset, Xvars, beta_hat, coeffCovar, colnums_bf, colnum_intercept)
% Function to compute the weighted combination of beta coefficients and its
% derivatives following a a FEMA analysis using spline basis functions
%
%% Inputs:
% basisSubset:      [b x q]     matrix of q basis functions output by
%                               createBasisFunctions.m
%
% Xvars:            [b x 1]     vector of X values over which basis
%                               functions were evaluated; output by
%                               createBasisFunctions.m
%
% beta_hat:         [p x v]     array of p beta estimates for v variables
%
% coeffCovar:       [p x p x v] estimated coefficient covariance matrix
%
% colnums_bf:       [numeric]   column number(s) of beta_hat corresponding
%                               to basis functions
%
% colnum_intercept: [numeric]   column number of beta_hat corresponding to
%                               the intercept
%
%% Outputs:
% beta_bf:          [b x v]     array of values showing fitted spline 
%                               trajectories for different values of Xvars
%                               (weighted sum of beta coefficients)
%
% SEmat_bf:         [b x v]     array of standard errors of mat_bf
%
% Zmat_bf:          [b x v]     array of Z scores of mat_bf
%
% beta_dbf:         [b x v]     array of values showing derivatives of
%                               spline estimates for different values of
%                               Xvars (weighted sum of betas for the
%                               derivatives of the splines)
%
% SEvols_dbf:       [b x v]     array of standard error of mat_dbf
%
% Zvols_dbf:        [b x v]     array of Z scores of mat_dbf
%
%% Notes:
% Adapted from FEMA_convert_splinevols but the ordering of both the inputs
% and the outputs is different

%% Initialize
num_X    = length(Xvars);
num_y    = size(beta_hat,2);

beta_bf  = NaN(num_X, num_y);
SE_bf    = NaN(num_X, num_y);
Z_bf     = NaN(num_X, num_y);

beta_dbf = NaN(num_X, num_y);
SE_dbf   = NaN(num_X, num_y);
Z_dbf    = NaN(num_X, num_y);

%% Calculate derivatives of basis functions
statmat_dbf = beta_hat(colnums_bf,:); % Select betas for columns of interest -- this shouldn't be needed
statmat_bf  = beta_hat([colnums_bf colnum_intercept],:);

[~, dbfmat] = gradient(basisSubset); 
dbfmat      = dbfmat/(Xvars(2)-Xvars(1)); % equivalent to gradient(basisSubset, unique(diff(Xvars)))

%% Compute spline function as weighted sum of betas for basis functions
% Gradient of basis functions
for xi = 1:num_X
    wvec            = dbfmat(xi,:)';
    valvec_dbf      = sum(statmat_dbf.*wvec,1); 
    sevec_dbf       = rowvec(sqrt(pagemtimes(pagemtimes(wvec',coeffCovar(colnums_bf,colnums_bf,:)),wvec)));
    zvec_dbf        = valvec_dbf ./ sevec_dbf;
    beta_dbf(xi, :) = valvec_dbf;
    SE_dbf(xi,:)    = sevec_dbf;
    Z_dbf(xi, :)    = zvec_dbf;
end

% Basis functions
for xi = 1:num_X
    wvec            = basisSubset(xi,:)';
    valvec_bf       = sum(statmat_bf.*[wvec; 1],1); 
    sevec_bf        = rowvec(sqrt(pagemtimes(pagemtimes([wvec; 1]',coeffCovar([colnums_bf colnum_intercept], [colnums_bf colnum_intercept],:)),[wvec; 1])));
    zvec_bf         = valvec_bf ./ sevec_bf;
    beta_bf(xi, :)  = valvec_bf;
    SE_bf(xi,:)     = sevec_bf;
    Z_bf(xi, :)     = zvec_bf;
end