function IDXout = IDX_H_std(IDX,data)
% IDXout = IDX_H_std(IDX,data) calculates the silhoutte values for each assignment pixel in IDX using data as a target.
disp('Calculating silhouette...')
if numel(size(data)) == 3
    sz = size(data);
    data = reshape(data,[prod(sz(1:2)) sz(3)]);
end
IDX(IDX==0) = NaN;
data(IDX==0,:) = NaN;
IDXout = reshape(silhouette(data, IDX(:),'correlation'),[256 256]);