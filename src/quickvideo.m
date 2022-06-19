function quickvideo(data,fr,cax,cmap,filename)
close all
figure('Position',[0 0 1000 500]);
v = VideoWriter(filename);
v.FrameRate = fr;
open(v)
for i = 1:size(data,3)
    imagesc(data(:,:,i)); axis image; axis off; colormap(cmap)
    colorbar
    caxis(cax)
    f = getframe(gca);
    writeVideo(v,f);
end
close(v);
close all