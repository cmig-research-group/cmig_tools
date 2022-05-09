function tfce_perm = FEMA_tfce(zmat_perm,colsinterest,datatype,mask)

logging('Running TFCE...');

nperms=size(zmat_perm,3);
tfce_perm = NaN(length(colsinterest),size(zmat_perm,2),size(zmat_perm,3));

for colsinteresti = 1:length(colsinterest)

      for permi = 1:nperms

            if mod(permi,100) == 0
                  logging('colsinteresti=%d permi=%d/%d',colsinteresti,permi,nperms);
            end

            zvec=zmat_perm(colsinteresti,:,permi);

            switch datatype
                  case 'vertex'
                        load('SurfView_surfs','icsurfs','surf_lh_pial','surf_rh_pial');
                        ic = NaN;
                        for ic = 1:length(icsurfs)
                              if 2*size(icsurfs{ic}.vertices,1) == length(zvec)
                                    icnum = ic;
                              end
                        end
                        icsurf = icsurfs{icnum};
                        nverts = size(icsurf.vertices,1);
                        nfacess = size(icsurf.faces,1);
                        statvec_lh = abs(zvec([1:nverts]));        data_lh = struct(); data_lh.fac = icsurf.faces; data_lh.vtx = surf_lh_pial.vertices(1:nverts,:);
                        statvec_rh = abs(zvec(nverts+[1:nverts])); data_rh = struct(); data_rh.fac = icsurf.faces; data_rh.vtx = surf_rh_pial.vertices(1:nverts,:);
                        statvec_tfce_lh = FEMA_tfce_vertex(statvec_lh,data_lh);
                        statvec_tfce_rh = FEMA_tfce_vertex(statvec_rh,data_rh);
                        statvec_tfce = [statvec_tfce_lh statvec_tfce_rh];
                  case 'voxel'
                        statvec = abs(zvec);
                        statvec_tfce = NaN(size(statvec));
                        statvol = fullvol(statvec,mask);
                        statvol_tfce = FEMA_tfce_voxel(statvol);
                        ivec_mask = find(mask>=0.5);
                        statvec_tfce(:) = statvol_tfce(ivec_mask);
                  otherwise
                        logging('%s: datatype=%s not supported',datatype);
                  return
            end

            tfce_perm(colsinteresti,:,permi)=statvec_tfce;

      end

end

end

