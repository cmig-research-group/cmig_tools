function residualVariable = FEMA_OLSResiduals(variable, fixedEffects, useFlag)
% OLS residualize a variable for other fixed effects
%% Inputs:
% variable:         [n x m]     matrix of n subjects and m varibles
% 
% fixedEffects:     [n x p]     matrix of n subjects and p covariates
% 
% useFlag:          character   useful for debugging; can specify which
%                               method to use for solving; can be:
%                                   * 'lsqminnorm'
%                                   * 'pinv'
%                                   * 'backslash'
%
%% Output(s):
% residualVariable  [n x m]     matrix of n subjects and m variables with 
%                               each variable residualized for the effect 
%                               of the fixedEffects using OLS
%
%% Notes:
% This is a simplified version of FEMA_residualizeGenotype without any
% involvement of bin - a variable (or a matrix of variables like genotype
% matrix) can be passed as inputs and a set of fixed effects to account for
% - the result is OLS residuals

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
end

%% Estimate beta
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

% Calculate beta: inv(X' * X) * X' * y
beta = iXtX * fixedEffects' * variable;

%% Compute residuals
residualVariable = variable - (fixedEffects * beta);