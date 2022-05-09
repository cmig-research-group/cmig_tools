function argout = ndgrid_amd(argin)

ndim = length(argin);
siz = NaN(1,ndim);
for i = 1:ndim
  siz(i) = numel(argin{i});
end

argout = NaN(prod(siz),ndim);
for i = 1:ndim
  x = argin{i}(:); % Extract and reshape as a vector.
  s = siz; s(i) = []; % Remove i-th dimension
  x = reshape(x(:,ones(1,prod(s))),[length(x) s]); % Expand x
  argout(:,i) = colvec(permute(x,[2:i 1 i+1:ndim])); % Permute to i'th dimension
end
