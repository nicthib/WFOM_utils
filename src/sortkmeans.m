function [IDXsorted,I] = sortkmeans(IDX)
c = []; inds = unique(IDX(:)); inds(inds == 0) = [];
inds(isnan(inds)) = [];
for i = 1:numel(inds)
    tmp = zeros(size(IDX));
    tmp(IDX==inds(i)) = 1; tmp = bwareaopen(tmp,100);
    if sum(tmp(:)) > 0
        stats = regionprops(tmp,'centroid');
        c(end+1,:) = stats.Centroid;
    else
        c(end+1,:) = [0 0];
    end
end
[~,I] = sort(c(:,2),'ascend');
IDXsorted = zeros(size(IDX));
for i = 1:size(I,1)
    IDXsorted(IDX==inds(I(i))) = i;
end