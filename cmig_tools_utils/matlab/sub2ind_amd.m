function ndx = sub2ind_amd(siz,submat)

%keyboard
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
for i = 1:length(siz)
    v = submat(:,i);
    ndx = ndx + (v-1)*k(i);
end
