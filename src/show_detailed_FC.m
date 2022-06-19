
function show_detailed_FC(FC,IDX,networks)
titles = 'ABCDEF';
if min(size(FC)) == 1
    FC = reshape(FC,[sqrt(numel(FC)),sqrt(numel(FC))]);
end
cols = numel(networks)+2;
networks = [networks networks];
regions = [];
for i = 1:numel(networks)
    addl = sum(networks(1:i-1));
    if isempty(addl); addl = 0; end
    regions{i} = [1:networks(i)]+addl;
end

subplot(2,cols,[1 2 1+cols 2+cols])
imagesc(FC); axis image; axis off; caxis([0 1])
%title(['State ' k])
add_group_lines(networks)
colorbar
for j = 1:numel(networks)
    spidx = j+2;
    if j > numel(networks)/2
        spidx = spidx+2;
    end
    subplot(2,cols,spidx)
    IDXtmp = assign_IDX(IDX,mean(FC(:,regions{j}),2));
    im = imagesc(IDXtmp); caxis([0 1])
    im.AlphaData = bwareaopen(IDX~=0,2000) & IDXtmp ~= 0;
    axis image; axis off; hold on
    highlight_network(IDX,regions{j},'k')
    
end
for j = 1:numel(networks)/2
    spidx = j+2;
    subplot(2,cols,spidx)
    title(titles(j))
end
colormap jet