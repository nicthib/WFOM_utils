function show_state_centroids(cc,row,col,networks)
% Regions
r = sqrt(size(cc,2));
cc = reshape(cc,[size(cc,1) r r]);
k = size(cc,1);
cc(cc==0) = NaN;
for j = 1:k
    subplot(row,col,j)
    imagesc(squeeze(cc(j,:,:))); hold on
    axis image; axis off;
    add_network_lines(networks)
    title(['S' mat2str(j)])
    colormap(jet)
    caxis([0 1])
end

% 
% figure
% title('Significant + correlations')
% for j = 1:k
%     axes(ha(j))
%     showsigcorrs(reshape(c(j,:),[92 92]),IDX,.9)
%     subplot(5,1,5)
%     hist(c(j,:),linspace(-1,1,201))
%     ylim([0 300])
%     axis off
% end
% 
% 
% figure
% %ha = tight_subplot(3,4);
% title('Correlation distribution')
% c(c==0) = NaN;
% for j = 1:k
%     subplot(3,4,j)
%     %axes(ha(j))
%    
% end
% 
% figure
% ha = tight_subplot(3,4);
% title('Significant + correlations')
% for j = 1:k
%     axes(ha(j))
%     showsigcorrs(reshape(-c(j,:),[92 92]),IDX,0)
% end