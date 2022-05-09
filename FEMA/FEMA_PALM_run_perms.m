function [zmat_perm beta_hat_perm beta_se_perm] = FEMA_PALM_run_perms(X,ymat,nperms,perms,contrast_vec,G_flag,colsinterest,VG)

      %Generates voxelwise permuted z statistics based on the specified
      %permutation matrix (perms)
      
      %INPUTS:
      %   X = design matrix (n obs by m predictors)
      %   ymat = brain matrix (n obs by v vertices)
      %   nperms = number of permutations (inc unpermuted as 1st col) - nperms is nperms given +1
      %   perms = permutation matrix (n obs by nperms) --> 1st column is unpermuted indices
      %   contrast_vec = m by 1 vector of
      
      %The following are generated from ABCD_PALM_produceperms.m --> ivec_perms, nperms, perms
      
      %OUTPUTS
      %   stat_ALL = mass univariate z stats across vertices for each permutation (nperms by v)
      
      %REFERENCE
      %Determining the association between cortical morphology and cognition in 10,145 children from the Adolescent Brain and Cognitive Development (ABCD) study using the MOSTest
      %C E Palmer, W Zhao, R Loughnan, C C Fan, W Thompson, T L Jernigan, A M Dale
      %bioRxiv 816025; doi: https://doi.org/10.1101/816025
      
      %(C) Written by Clare E Palmer & Anders Dale, Sept 2019
      
      if nperms<size(perms,2)
            nperms=nperms+1;
      end
      
      if isempty(contrast_vec)
            contrast_all=eye(size(X,2));
            contrast_all=contrast_all(colsinterest,:)';
            beta_hat_perm = NaN([size(X,2) size(ymat,2) nperms]);
            beta_se_perm = NaN([size(X,2) size(ymat,2) nperms]);
            zmat_perm = NaN([size(X,2) size(ymat,2) nperms]);
      elseif ~isempty(contrast_vec)
            colsinterest=[1,colsinterest+1];
            contrast_all=nan(size(X,2)+1,size(X,2));
            contrast_all(1,:)=contrast_vec;
            contrast_all(2:end,:)=eye(size(X,2));
            contrast_all=contrast_all(colsinterest,:)';
            beta_hat_perm = NaN([size(X,2)+1 size(ymat,2) nperms]);
            beta_se_perm = NaN([size(X,2)+1 size(ymat,2) nperms]);
            zmat_perm = NaN([size(X,2)+1 size(ymat,2) nperms]);;
      end
      
      contrast_vec=colvec(contrast_vec);

      M=X;

      for cols=1:size(contrast_all,2)
                  
            Xpart=X*contrast_all(:,cols);
            Zpart=X(:,(contrast_all(:,cols)==0));

            nvols = size(M,1);
            nparams = size(M,2);
            
            Mi = pinv(M);
            MiM = Mi*M; % Crosstalk matrix
            mim = diag(MiM);
            Mi(find(mim<0.9),:) = 0; % Eliminate non-identifiable params
            C = Mi*Mi';
            
            I=eye(size(Zpart,1));
            ez=(I-(Zpart*pinv(Zpart)))*ymat;
            
            fprintf(1,'START PERMS (colsinterest, %d): 1 unpermuted and %d perms (now=%s)\n',colsinterest(cols),nperms-1,datestr(now));
            start_t=now();
            for iter = 1:nperms

                  if mod(iter,10)==0
                        fprintf(1,'iter=%d of %d (now=%s)\n',iter,nperms,datestr(now));
                  end
            
                  ivec_tmp=perms(:,iter);
            
                  psi=Mi*ez(ivec_tmp,:);
                  res=(I-M*Mi)*ez(ivec_tmp,:);

                  if G_flag==0

                        vol_sum2 = sum(res.^2,1);
                        sig_vol = sqrt(vol_sum2/(nvols-nparams));

                        zvec=(psi'*contrast_all(:,cols))'./(sqrt(C(1,1))*sig_vol);
                        zvec(~isfinite(zvec)) = NaN; %Need to mask imaging data
                        se=sqrt(C(1,1))*sig_vol;

                        zmat_perm(colsinterest(cols),:,iter)=zvec; %predictor x vertex/voxel x perm
                        beta_se_perm(colsinterest(cols),:,iter)=se; %predictor x vertex/voxel x perm
                        beta_hat_perm(colsinterest(cols),:,iter)=psi'*contrast_all(:,cols); %predictor x vertex/voxel x perm


                  elseif G_flag==1

                        if iter==1
                              fprintf(1,'Calculating winkler pivotal, G stat.  Will take long time. (now=%s)\n',datestr(now));
                              if isempty(VG)
                                    fprintf(1,'No VG provided. Using default. (now=%s)\n',datestr(now));
                              end
                        end

                        %Due to speed only runs G for contrast_vec
                        G = winkler_pivotal(M,psi,res,contrast_all(:,cols),VG,ymat);
                        zmat_perm(colsinterest(cols),:,iter)=G;
                        beta_hat_perm(colsinterest(cols),:,iter)=psi'*contrast_all(:,cols);

                        fprintf(1,'iter=%d of %d (now=%s)\n',iter,nperms,datestr(now));

                        beta_se_perm=[];

                  end
            
            end
            end_t=now();
            fprintf(1,'END PERMS: Time taken to complete 1 unpermuted and %d perms = %0.2f seconds (now=%s)\n',nperms-1,(end_t-start_t)*3600*24,datestr(now));

      end
      
      
end
      
      
      