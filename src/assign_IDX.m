function IDXout = assign_IDX(IDX,val)
IDXout = zeros(size(IDX));
for i = 1:numel(val)
    IDXout(IDX==i) = val(i);
end