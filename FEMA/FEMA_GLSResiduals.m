function [residualVariable, XtW, Bi, useFlag, useLSQ, lowRank] = FEMA_GLSResiduals(variable, fixedEffects, XtW, Bi, allWsFam, useFlag, useLSQ, lowRank)
% GLS residualize a variable for other fixed effects
%% Inputs:
% variable:         [n x m]     matrix of n subjects and m variables
% 
% fixedEffects:     [n x p]     matrix of n subjects and p covariates
%
% XtW:              [p x p]     matrix containing the inverse of the
%                               X' * inv(V) * X term
%
% Bi:               [p x p]     inverse of the X' * inv(V) * X term
%
% allWsFam:         [n x n]     sparse matrix containing compiled inverse
%                               of the V term (see, FEMA_compileTerms)
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
% XtW:              [p x n]     matrix containing the X' * inv(V) term,
%                               which can be passed back for faster
%                               computation (for example, when
%                               residualizing multiple SNPs fopr the same X
%                               variables)
% 
% Bi:               [p x p]     matrix containing the inv(X' * inv(V) * X) 
%                               term, which can be passed back for faster
%                               computation (for example, when
%                               residualizing multiple SNPs fopr the same X
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
% as inputs and a set of fixed effects to account for, along with the
% variance-covariance random effects matrix - GLS residuals are returned
%
% Either allWsFam or XtW and Bi need to be passed as inputs

%% Check if debugging mode
if ~exist('useFlag', 'var') || isempty(useFlag)
    % Determine if lsqminnorm can be used
    if ~exist('useLSQ', 'var') || isempty(useLSQ)
        if exist('lsqminnorm', 'file')
            useLSQ = true;
        else
            useLSQ = false;
        end
    end

    % Check if rank deficit
    if ~exist('lowRank', 'var') || isempty(lowRank)
        if rank(fixedEffects) < size(fixedEffects, 2)
            lowRank = true;
        else
            lowRank = false;
        end
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

%% Check if variables can be reused
if (exist('XtW', 'var') && ~isempty(XtW)) && ...
   (exist('Bi',  'var') && ~isempty(Bi))
    toReuse = true;
else
    toReuse = false;
end

%% Estimate beta coefficient
% Reuse, if required variables exist
if ~toReuse
    if ~exist('allWsFam', 'var') || isempty(allWsFam)
        error('Either provide allWsFam or XtW and Bi as inputs');
    end

    % Compute XtW
    XtW = double(fixedEffects)' * allWsFam;

    % Compute XtWX
    B  = XtW * fixedEffects;

    % Calculate inverse
    switch useFlag
        case 'lsqminnorm'
            Bi = lsqminnorm(B, eye(size(B)));
    
        case 'pinv'
            Bi = pinv(B);
    
        case 'backslash'
            Bi = B \ eye(size(B));
    end
end

% Calculate beta: inv(X' * inv(V) * X) * X' * inv(V) * y
beta = Bi * (XtW * variable);

%% Compute residuals
residualVariable = variable - (fixedEffects * beta);