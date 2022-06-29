function show_detailed_FCs(FCs,IDX,networks,final_ord)
% Iterated version of show_detailed_FC(), which shows all states in matrix FCs (nx92x92) in one figure.
titles = 'ABCDEF';
cols = numel(networks)+2;
k = size(FCs,1);
networks = [networks networks];
regions = [];
for i = 1:numel(networks)
    addl = sum(networks(1:i-1));
    if isempty(addl); addl = 0; end
    regions{i} = [1:networks(i)]+addl;
end

for F=1:k
    FC = reshape(FCs(F,:),[sqrt(size(FCs,2)),sqrt(size(FCs,2))]);
    FC = FC(final_ord,final_ord);
    spidx = 1+(F-1)*cols*2;
    subplot(2*k,cols,[spidx spidx+1 spidx+cols spidx+cols+1])
    imagesc(FC); axis image; axis off; caxis([0 1])
    title(['State ' num2str(F)])
    colorbar
    
    for j = 1:numel(networks)
        spidx = 2+(F-1)*cols*2+j;
        if j > numel(networks)/2
            spidx = spidx+2;
        end
        subplot(2*k,cols,spidx)
        IDXtmp = assign_IDX(IDX,mean(FC(:,regions{j}),2));
        im = imagesc(IDXtmp); caxis([0 1])
        im.AlphaData = bwareaopen(IDX~=0,2000) & IDXtmp ~= 0;
        axis image; axis off; hold on
        highlight_network(IDX,regions{j},'k')
    end
%     for j = 1:numel(networks)/2
%         spidx = j+2;
%         subplot(2,cols,spidx)
%         title(titles(j))
%     end
end
colormap jet