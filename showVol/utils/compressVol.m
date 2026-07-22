function c = compressVol(vol)
	%% compression/expansion
	%simple compression method: only store non-zero values, as single precision. since ROIs are very sparse
	%for the most part, the space savings (thus load time improvement) is
	%large. E.g. aseg only 1.4% of voxels are non-zero
	s = size(vol);
	if prod(s) >= 2^32, error('Volume too large to compress with uint32 indexes'), end
	nz = cast(find(vol(:)~=0),'uint32'); %index of non-zero voxels

	c = {s nz(:) single(vol(nz))}; %volume size, non-zero-voxel indices, non-zero-voxel values

	fprintf('(compressed to %s%% of original)\n',num2str(100*2*length(nz)/prod(s),2)); %factor of 2 since we save index and values

	% companion expansion function when loading atlas. Not used here, but FYI.
	%  e.g. 
	% if iscell(aseg.prob)
	%     aseg.prob = expandVol(aseg.prob);
	% end
end 