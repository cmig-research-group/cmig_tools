function [Ws_famtype, Ws_fam] = FEMA_compileTerms(clusterinfo, binvec, nfamtypes, famtypevec, sig2mat, RandomEffects, GroupByFamType, SingleOrDouble, OLSflag)
%% Compile Vs and Ws terms for every bin
%% Inputs:
% clusterinfo:      [1 x f]     cell type of f families/clusters, returned
%                               as an output from FEMA_parse_family
% 
% binvec:           [1 x b]     vector of b bin values
% 
% nfamtypes:        [1 x 1]     used if GroupByFamType
% 
% famtypevec:       [1 x f]     vector of family types
% 
% sig2mat:          [r x b]     estimates of r random effects across b bins
% 
% RandomEffects:    [1 x r]     cell having names of the r random effects
% 
% GroupByFamType:   logical     whether to GroupByFamType or not
% 
% SingleOrDouble:   character   single or double precision
% 
% OLSflag:          logical     OLS or GLS solution
%
%% Outputs:
% Ws_famtype and Ws_fam as cell types containing inverses of the "V" term 
% across clusters and bins

%% Initialize
nfam         = length(clusterinfo);
nbins        = length(unique(binvec(isfinite(binvec)),'stable'));
Vs_famtype   = cell(nbins,nfamtypes); 
Vis_famtype  = cell(nbins,nfamtypes); 
Ws_famtype   = cell(nbins,nfamtypes);
Vs_fam       = cell(nbins,nfam); 
Vis_fam      = cell(nbins,nfam); 
Ws_fam       = cell(nbins,nfam);
binLoc       = 1;

for sig2bini = unique(binvec(isfinite(binvec)),'stable')
    ivec_bin = find(binvec==sig2bini);
    sig2vec  = mean(sig2mat(:,ivec_bin), 2);
    
    if ~isempty(ivec_bin)
        if GroupByFamType % Compute Vs and Vis by family type
            for fi = 1:nfamtypes
                ivec = find(famtypevec==fi);
                Vs_famtype{binLoc, fi} = 0;
                for ri = 1:length(RandomEffects)
                    Vs_famtype{binLoc, fi} = Vs_famtype{binLoc, fi} + sig2vec(ri)*clusterinfo{ivec(1)}.(['V_', RandomEffects{ri}]);
                end
                Vis_famtype{binLoc, fi} = cast(double(Vs_famtype{binLoc, fi}) \ eye(size(Vs_famtype{binLoc, fi})), SingleOrDouble);
                % Vis_famtype{binLoc, fi} = cast(pinv(double(Vs_famtype{binLoc, fi})),SingleOrDouble);
                if OLSflag
                    Ws_famtype{binLoc, fi} = eye(size(Vis_famtype{binLoc, fi}), class(Vis_famtype{binLoc, fi}));
                else
                    Ws_famtype{binLoc, fi} = Vis_famtype{binLoc, fi};
                end
            end
        else % Compute Vs and Vis for each family
            for fi = 1:nfam
                Vs_fam{binLoc, fi} = 0;
                for ri = 1:length(RandomEffects)
                    Vs_fam{binLoc, fi} = Vs_fam{binLoc, fi} + sig2vec(ri)*clusterinfo{fi}.(['V_', RandomEffects{ri}]);
                end
                Vis_fam{binLoc, fi} = cast(double(Vs_fam{binLoc, fi}) \ eye(size(Vs_fam{binLoc, fi})), SingleOrDouble);
                % Vis_fam{binLoc, fi} = cast(pinv(double(Vs_fam{binLoc, fi})),SingleOrDouble);
                if OLSflag
                    Ws_fam{binLoc, fi} = eye(size(Vis_fam{binLoc, fi}), class(Vis_fam{binLoc, fi}));
                else
                    Ws_fam{binLoc, fi} = Vis_fam{binLoc, fi};
                end
            end
        end
    end
    binLoc = binLoc + 1;
end