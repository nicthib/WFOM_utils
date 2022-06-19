function IDX_n = IDXtodigraph(IDX)
IDX_n = [];
for i = 1:max(IDX(:))
    tmp = IDX==i;
    tmp = imdilate(tmp,ones(3));
    neighbors = unique(IDX(find(tmp)));
    neighbors(isnan(neighbors)) = [];
    neighbors(neighbors==i) = [];
    IDX_n{i} = neighbors;
end