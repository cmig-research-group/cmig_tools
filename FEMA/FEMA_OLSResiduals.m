function [residualVariable, iXtX, useFlag, useLSQ, lowRank] = FEMA_OLSResiduals(variable, fixedEffects, iXtX, useFlag, useLSQ, lowRank)
% OLS residualize a variable for other fixed effects
%% Inputs:
% variable:         [n x m]     matrix of n subjects and m variables
% 
% fixedEffects:     [n x p]     matrix of n subjects and p covariates
% 
% iXtX:             [m x m]     matrix containing the inverse of the X'X
%                               term (optional)
%
% useFlag:          character   can specify which method to use for
%                               computing inverse; can be: 
%                                   * 'lsqminnorm'
%                                   * 'pinv'
%                                   * 'backslash'
%
% useLSQ:           logical     true or false specifying if lsqminnorm can
%                               be used
%
% lowRank:          logical     true or false specifying if the
%                               fixedEffects matrix is rank deficient
% 
%% Output(s):
% residualVariable: [n x m]     matrix of n subjects and m variables with 
%                               each variable residualized for the effect 
%                               of the fixedEffects using OLS
%
% iXtX:             [p x p]     matrix containing the inverse of the X'X
%                               term, which can be passed back for faster
%                               computation (for example, when
%                               residualizing multiple SNPs for the same X
%                               variables)
%
% useFlag:          character   the method used for computing inverse
% 
% useLSQ:           logical     whether lsqminnorm can be used
% 
% lowRank:          logical     whether fixedEffects matrix is rank deficit
%
%% Notes:
% A variable (or a matrix of variables like genotype matrix) can be passed
% as inputs and a set of fixed effects to account for - OLS residuals are
% returned

%% Check if debugging mode
if ~exist('useFlag', 'var') || isempty(useFlag)    
    % Determine if lsqminnorm can be used
    if exist('lsqminnorm', 'file')
        useLSQ = true;
    else
        useLSQ = false;
    end

    % Check if rank deficit
    if rank(fixedEffects) < size(fixedEffects, 2)
        lowRank = true;
    else
        lowRank = false;
    end

    if lowRank
        if useLSQ
            useFlag = 'lsqminnorm';
        else
            useFlag = 'pinv';
        end
    else
        if useLSQ
            useFlag = 'lsqminnorm';
        else
            useFlag = 'backslash';
        end
    end
else
    useFlag = lower(useFlag);
    if ~ismember(useFlag, {'lsqminnorm', 'backslash', 'pinv'})
        error(['Unknown useFlag provided: ', useFlag, '; should be one of: lsqminnorm, backslash, or pinv']);
    end
end

%% Estimate beta
% Reuse iXtX, if it exists
if ~exist('iXtX', 'var') || isempty(iXtX)
    XtX  = fixedEffects' * fixedEffects;

    % Calculate inverse
    switch useFlag
        case 'lsqminnorm'
            iXtX = lsqminnorm(XtX, eye(size(XtX)));
    
        case 'pinv'
            iXtX = pinv(XtX);
    
        case 'backslash'
            iXtX = XtX \ eye(size(XtX));
    end
end

% Calculate beta: inv(X' * X) * X' * y
beta = iXtX * (fixedEffects' * variable);

%% Compute residuals
residualVariable = variable - (fixedEffects * beta);