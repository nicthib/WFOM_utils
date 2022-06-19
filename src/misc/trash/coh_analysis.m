close all
for i = 1:50
subplot(121)
imagesc(reshape(cv(:,i),[36 36]))
caxis([0 1])
axis image
colormap jet
colorbar
title('Peak coherence value')

subplot(122)
imagesc(reshape(cf(:,i),[36 36]))
caxis([0 10])
axis image
colorbar
title('Peak cohesion frequency')
pause(.1)
end