function [ord,IDX_sub,FCIDX_sort] = FC_cluster(H,IDX,k)
show = 1;
% sum bilaterally
H1 = H(1:size(H,1)/2,:);
H2 = H(size(H,1)/2+1:size(H,1),:);
H_nb = H1+H2;
FC = corr(H_nb');
FCIDX = kmeans(FC,k,'Replicates',100); % FC ordering

IDX_sub = zeros(size(IDX));
for i = 1:k
   tmp = ismember(IDX,find(FCIDX==i));
   IDX_sub(tmp==1) = i;
end
[IDX_sub, subregion_ord] = sortkmeans(IDX_sub);

FCIDX_sort = zeros(size(FCIDX));
for i = 1:k
    FCIDX_sort(FCIDX == subregion_ord(i)) = i;
end
[~, ord] = sort(FCIDX_sort);

if show
    figure
    nanrow = nan(numel(FCIDX_sort),1);
    subplot(131)
    i1 = imagesc([FCIDX_sort/k nanrow FC]);
    axis image; axis off; colormap jet; caxis([0 1])
    title('Input FC map')
    
    subplot(132);
    i2 = imagesc([FCIDX_sort(ord)/k nanrow FC(ord,ord)]);
    axis image; axis off; caxis([0 1])
    title('Sorted FC map')

    subplot(133); imagesc(IDX_sub)
    axis image; axis off;
    title(['Brain regions for k = ' mat2str(k)])
    
    i1.AlphaData = ~isnan(i1.CData);
    i2.AlphaData = ~isnan(i2.CData);
end

ord = [ord; ord+numel(ord)];
FCIDX_sort = [FCIDX_sort;FCIDX_sort];