function v = rowvec(M)

% v = reshape(M,[1 prod(size(M))]);
v = reshape(M,[1 numel(M)]);