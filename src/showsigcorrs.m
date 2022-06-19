function showsigcorrs(c,IDX,thr,color,n)
c = reshape(c,[sqrt(numel(c)) sqrt(numel(c))]);
c = tril(c);
c(c>1) = 1;
c(c<thr) = 0;
r = sqrt(numel(c));
[~,I] = sort(c(:),'ascend');
[sigcx,sigcy] = ind2sub([r,r],I);
IDX(isnan(IDX)) = 0;
%IDXg = double(imerode(logical(imgradient_nic(IDX)),ones(2))); %axis image; axis off; colormap gray; hold on

%showIDX_labels_subregions(IDX,[9 16 9 5 7]); hold on
%im = imagesc(IDXg); axis image; axis off; hold on; caxis([0 2])
%im.AlphaData = IDX~= 0;
nlines = 0; i=1;
while nlines < n 
    if c(I(i)) > thr
        tmp1cc = bwconncomp(bwareaopen(IDX==sigcx(i),50));
        tmp2cc = bwconncomp(bwareaopen(IDX==sigcy(i),50));
        for j = 1:numel(tmp1cc.PixelIdxList)
            for k = 1:numel(tmp2cc.PixelIdxList)
                tmp1 = zeros(size(IDX));
                tmp2 = zeros(size(IDX));
                tmp1(tmp1cc.PixelIdxList{j}) = 1;
                tmp2(tmp2cc.PixelIdxList{k}) = 1;
                s1 = regionprops(tmp1); c1 = s1.Centroid;
                s2 = regionprops(tmp2); c2 = s2.Centroid;
                dist = sqrt((c2(1)-c1(1))^2+(c2(2)-c1(2))^2);
                %if dist > 50
                    ctmp = c(I(i));
                    %color = round(interp1(cax,[0 100],ctmp));
                    line([c1(1) c2(1)],[c1(2) c2(2)],'LineWidth',1,'Color',color)
                    nlines = nlines+1;
                %end
            end
        end
    end
    i=i+1;
    if i > numel(c)
        break
    end
end
