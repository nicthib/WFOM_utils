% cleanupIDX(IDX,a) removes small noisy components from a ROI map. a defines the size (in pixels) of spurious components to remove.
function out = cleanupIDX(IDX,a)
out = zeros(size(IDX));
for i = 1:max(IDX(:))
    tmp = bwareaopen(IDX == i,a);
    tmp = tmp*i;
    out = out+double(tmp);
end
    