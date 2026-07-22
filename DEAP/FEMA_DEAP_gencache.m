function FEMA_DEAP_gencache(nfrac, dirname_cache, dirname_imaging, dirname_out, fname_design)
% fname_design
% dirname_voxelwise - dependant on the data release v4, v5, v6, etc...
% nfrac - need to pass to init
%
% dirname_out we can remove?
%

if isdeployed
    nfrac = str2double(nfrac);
    if isnan(nfrac)
        error('One or more input arguments could not be converted to a number.');
    end
end

logging('FEMA_DEAP_gencache Begin %d, %s, %s, %s, %s', nfrac, dirname_cache, ...
        dirname_imaging, dirname_out, fname_design);

% Print build/compile information
% print_compile_stats();
% Cache voxelwise phenotypes
% This need not be hardcoded at all - should be able to meaningfully parse
% the ABCD release directories and subdirectories
%  measlist = {{'rsi' 'RNI'} {'rsi' 'RND'} {'smri' 'JA'}}; % List of phenotypes to cache
% measlist = {{'smri' 'JA'}};
measlist = {{'smri' 'thk_sm16'}};
datatype = 'vertex';
for measi = 1:length(measlist)
    subdir = measlist{measi}{1};
    voxelwise_pheno = measlist{measi}{2};
    fstem_pheno_voxel = sprintf('%s', voxelwise_pheno);
    if exist(dirname_out, 'dir')
        warning('Directory already exists; overwriting');
    else
        % jsystem(sprintf('rm -rf %s',dirname_out));
        mkdir(dirname_out);
    end
    dirname_voxelwise_tmp = sprintf('%s/%s', dirname_imaging, subdir);

    % Should add dirname_cache and nfrac as arguments, get rid of dirname_out? 
    % -- should write out imaging phenotypes only
    FEMA_DEAP_init(fstem_pheno_voxel, fname_design, dirname_out, dirname_voxelwise_tmp, ...
                   datatype, dirname_cache, nfrac);
    logging('***Done*** with %s', voxelwise_pheno);
end
logging('FEMA_DEAP_gencache Complete');
end
