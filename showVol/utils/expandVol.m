function vol = expandVol(c)\
	%% compression/expansion
	%simple compression method: only store non-zero values, as single precision. since ROIs are very sparse
	%for the most part, the space savings (thus load time improvement) is
	%large. E.g. aseg only 1.4% of voxels are non-zero
	
	vol = zeros(c{1},'single');
	vol(c{2}) = c{3};
end 