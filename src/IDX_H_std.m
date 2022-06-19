function IDXout = IDX_H_std(IDX,data)
disp('Calculating silhouette...')
if numel(size(data)) == 3
    sz = size(data);
    data = reshape(data,[prod(sz(1:2)) sz(3)]);
end
%[C,S] = fsvd(data,500,2,1);
IDX(IDX==0) = NaN;
data(IDX==0,:) = NaN;
IDXout = reshape(silhouette(data, IDX(:),'correlation'),[256 256]);
