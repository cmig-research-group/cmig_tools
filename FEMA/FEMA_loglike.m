function nlogl = FEMA_loglike(x,yvec_res,clusterinfo,RandomEffect_list)

disp(x)

nlogl = 0;
for fi = 1:length(clusterinfo)
  Sigma = 0;
  for ri = 1:length(RandomEffect_list)
    Sigma = Sigma + x(ri)*getfield(clusterinfo{fi},sprintf('V_%s',RandomEffect_list{ri}));
  end
  nlogl = nlogl + -log(mvnpdf(double(yvec_res(clusterinfo{fi}.jvec_fam)),0,double(Sigma)));
end
%disp(x)
%disp(nlogl)
