function tfcestat = FEMA_tfce_vertex(D,data)

Yarea= palm_calcarea(data,1);
Yadjacency = palm_adjacency(data.fac,1);

% Define default "delta h". Can also be specified by the user
dh = max(D(:))/100;
E = 0.5;
H = 2.0;

%Compute the TFCE statistic
tfcestat = zeros(size(D));
for h = dh:dh:max(D(:))
  %Label a DPV file using adjacency matrix
  dpxl  = palm_dpxlabel(D>=h,Yadjacency);
  U     = unique(dpxl(dpxl>0))';
  for u = 1:numel(U)
    idx = dpxl == U(u);
    tfcestat(idx) = tfcestat(idx) + sum(Yarea(idx)).^E * h^H;
  end
end

%Correct for delta h
tfcestat = tfcestat * dh;

