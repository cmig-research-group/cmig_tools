
function [delta_r2] = FEMA_calcdeltaR2(FEMA_outfile_reduced,FEMA_outfile_full, ymat)

% Delta R2 = difference between variance explained of full compared to reduced model
      
% Calculates average voxelwise delta R2 between full and reduced model using SSE output from a full and reduced model
% Average numerator and denominator to avoid errors dividing by zero that may occur at some voxels

      % y = imaging data used to compute model (ymat)
      % y_hat = predicted imaging data from model
      % mean R2_full = mean(var(y_hat_full)) / mean(var(y))
      % mean R2_red = mean(var(y_hat_red)) / mean(var(y))

      % Mean voxelwise delta R2 = R2_full - R2_red

% Requires ymat as input - total variance is equal to var(ymat)
      
% INPUTS:
%     FEMA_outfile_full = path to SSE_wrapper or SSE_fit output for FULL model
%     FEMA_outfile_reduced = path to SSE_wrapper or SSE_fit output for REDUCED model
%     ymat = matrix of voxelwise imaging data used for modelling

      for model=1:2

            if model==1
                  load(FEMA_outfile_full);
            elseif model==2
                  load(FEMA_outfile_reduced);
            end
      
            if size(ymat,1)~=size(X,1)
      
                  error('Ymat and model output have different number of subjects. ymat should be the same data used to generate the models. The full and reduced models should be generated using the same data.')
      
            else

                  if size(colnames_model,2)>size(X,2)
                        col_diff=size(colnames_model,2)-size(X,2);
                  end
      
                  if exist('beta_hat')

                        betas=beta_hat;
                        betas(1:col_diff,:)=[];
                        betas=betas';

                  else

                        mask=vol_sig2bin;
                        mask(find(mask==0))=nan;
                        maskvec=reshape(mask, 100*100*130, 1);
      
                        betas=reshape(vol_beta_hat, 100*100*130,size(vol_beta_hat,4));
                        betas(isnan(maskvec),:)=[];
                        betas(:,1:col_diff)=[];

                  end

                  yhat=betas*X';  %beta = nvoxel x nscanner,  Xscan = nsubj * nscanner (dummy coded),  yhat = nvoxel x nsubj
      
                  sigma_fixed=var(yhat,[],2); %expected yhat dimensions nvoxels x nsubj
                  voxmean_fixed=mean(sigma_fixed);
      
                  sigma_y=var(ymat); %expected ymat dimensions nsubj x nvoxels
                  voxmean_dat=mean(sigma_y);
      
                  if model==1
                        full_r2=voxmean_fixed / voxmean_dat;
                  elseif model==2
                        red_r2=voxmean_fixed / voxmean_dat;
                  end

      
            end
      
      end

delta_r2 = full_r2 - red_r2;
      
end
      