function out = cleanupIDX(IDX,a)
out = zeros(size(IDX));
for i = 1:max(IDX(:))
    tmp = bwareaopen(IDX == i,a);
    tmp = tmp*i;
    out = out+double(tmp);
end
    