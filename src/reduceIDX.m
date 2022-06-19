function IDXout = reduceIDX(IDX)
inds = unique(IDX(:));
IDXout = zeros(size(IDX));
for i = 1:numel(inds)
    IDXout(IDX==inds(i)) = i;
end
IDXout = IDXout - 1;