function img = expandFOD(a)
if isempty(a{3})
    img = zeros(a{1},a{2},3,'uint8');
else
    img = zeros(a{1}*a{2},3,'uint8');
    img(a{3},:)= a{4};
    img = reshape(img,a{1},a{2},3);
end