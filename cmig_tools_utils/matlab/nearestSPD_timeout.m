function [Ahat, converge] = nearestSPD_timeout(A, maxIter)
% nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
% usage: Ahat = nearestSPD(A)
%
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% arguments: (input)
%  A - square matrix, which will be converted to the nearest Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.
%
% This modified version has a max number of iterations; if for any reason
% the nearest SPD has not been found by then, return the version that was
% available at that iteration and set converge to false - prevents infinite
% loop - PP (24 Feb 2025)

% if nargin ~= 1
%   error('Exactly one argument must be provided.')
% end

% Max iter is set to 100
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 100;
end

% test for a square matrix A
[r,c] = size(A);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
  % A was scalar and non-positive, so just return eps
  Ahat = eps;
  return
end

% symmetrize A into B
B = (A + A')/2;

% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[U,Sigma,V] = svd(B);
H = V*Sigma*V';

% get Ahat in the above formula
Ahat = (B+H)/2;

% ensure symmetry
Ahat = (Ahat + Ahat')/2;

% test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
p = 1;
k = 0;
iter = 0;
converge = true;
while p ~= 0
  iter = iter + 1;
  [R,p] = chol(Ahat);
  k = k + 1;
  if p ~= 0
    % Ahat failed the chol test. It must have been just a hair off,
    % due to floating point trash, so it is simplest now just to
    % tweak by adding a tiny multiple of an identity matrix.
    mineig = min(eig(Ahat));
    Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
  end
  if iter >= maxIter
      converge = false;
      break
  end
end