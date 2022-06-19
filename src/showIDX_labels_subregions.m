function showIDX_labels_subregions(IDX,subregions,cmap,bilat,ax)
nc = max(IDX(:));
if ~bilat; nc=nc*2; end
cmap = [0 0 0;cmap];
cIDX = zeros(256,256);
max_comp = 0;
for j = 1:numel(subregions)
    target_comps = max_comp+1:max_comp+subregions(j);
    cIDX(ismember(IDX,[target_comps target_comps+nc/2])) = j;
    max_comp = max_comp+subregions(j);
end

% show map
imagesc(cIDX);
axis image; axis off; colormap(ax,cmap);
im = getframe(gca); im = imresize(im.cdata,[256 256]);

drawnow; pause(.5); drawnow
im = imoverlay(im,imerode(logical(imgradient_nic(IDX)),ones(2)),'black');
% Show Original
im = imagesc(im); axis image; axis off;
im.AlphaData = IDX > 0;
for j = 1:max(IDX(:))
    tmp = IDX==j; tmp = bwareaopen(tmp,50);
    c = regionprops(tmp,'Centroid');
    for k = 1:numel(c)
        if j <= nc/2 || ~bilat
            txt = mat2str(j);
        else
            txt = mat2str(j-nc/2);
        end
        text(c(k).Centroid(1),c(k).Centroid(2),txt,'Color','w','HorizontalAlignment','center','VerticalAlignment','middle')
    end
end

