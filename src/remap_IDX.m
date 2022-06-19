function IDXout = remap_IDX(IDX,ord)
IDXout = zeros(size(IDX));
for i = 1:max(IDX(:))
    IDXout(IDX==ord(i)) = i;
end