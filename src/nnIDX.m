function IDX = nnIDX(IDX) % Nearest neighbor IDX
%IDXf = imfill(IDX,'holes');
%IDX(IDXf==0) = NaN;
[zr,zc] = find(IDX == 0);
[nzr, nzc] = find(IDX > 0);
for i = 1:numel(zr)
    [~,idx] = min(sum(([nzr,nzc]-repmat([zr(i),zc(i)],[numel(nzr),1])).^2,2));
    IDX(zr(i),zc(i)) = IDX(nzr(idx),nzc(idx));
end
