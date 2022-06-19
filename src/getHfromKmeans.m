% Takes an indexed map (IDX) and extracts timecourses from 3D datasets that
% correspond to those regions. Also erodes if desired
function [H,Hstd] = getHfromKmeans2(data,IDX,n_erode)
ss = size(data);
if numel(ss) == 3
    data = reshape(data,[prod(ss(1:2)) ss(3)]);
end
for i = 1:max(IDX(:))
    tmp = imfill(IDX == i,'holes');
    %tmp = bwareaopen(tmp,50);
    if n_erode > 0
        for j = 1:n_erode
            tmp = imerode(tmp,ones(3));
        end
    end
    tmp = double(tmp);
    tmp(tmp == 0) = NaN;
    tmp = tmp(:);
    H(i,:) = nanmean(data(~isnan(tmp),:));
    Hstd(i,:) = nanstd(data(~isnan(tmp),:));
end