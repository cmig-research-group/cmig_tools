function v = colvec(M)

% v = reshape(M,[prod(size(M)) 1]);
v = reshape(M,[numel(M) 1]);