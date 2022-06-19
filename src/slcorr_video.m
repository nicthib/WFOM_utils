function slcorr_video(H,st,c,cc,win,IDX,vid,savename)
ww = numel(win);
r = size(H,1);
k = size(c,1);
vO = VideoWriter(savename);
vO.FrameRate = 20;
open(vO);
vI = VideoReader(vid);
close all
figure('Position',[100 100 1200 600]); colormap jet
divs=[46 50];
for i = 1:1200
    subplot(5,4,1:8); cla; title(savename)
    imagesc(H(:,1:1200)); caxis([-.1 .1]); colorbar;
    set(gca,'XTickLabel',[]);
    line([i i],[0 size(H,1)+1],'Color','k','LineWidth',2)
    line([i+ww i+ww],[0 size(H,1)+1],'Color','k','LineWidth',2)
    ylabel('Component #')
    subplot(5,4,[9 9+4]) % Raw corr
    tmp = divcorrmap(reshape(cc(:,i),[r r]),divs);
    imagesc(tmp);
    im = imagesc(tmp); caxis([0 1]); axis off; axis image;
    im.AlphaData = ~isnan(tmp);
    title(sprintf('Time: %0.2f sec',i/20))
    subplot(5,4,[10 10+4]) % State label
    tmp = divcorrmap(reshape(c(st(i),:),[r r]),divs);
    im = imagesc(tmp); im.AlphaData = ~isnan(tmp);
    caxis([0 1]); title(st); axis off; axis image
    title(['State # ' mat2str(st(i))])
    subplot(5,4,[12 12+4])
    for j = 1:3
        mfr = readFrame(vI);
    end
    imagesc(mfr); axis image; axis off
    subplot(5,4,[11 11+4])
    showsigcorrs(reshape(cc(:,i),[r r]),IDX,.8)
    title('Top 100 connections > .8')
    
    subplot(5,4,17:20); cla
    plot(st(1:i),'Color','k','LineWidth',2); ylim([0 k+1]); xlim([0 1200])
    set(gca, 'box','off','XTickLabel',[],'XTick',[],'XColor','none','YTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12'},'YTick',1:k)
    currFrame = getframe(gcf);
    writeVideo(vO,currFrame);
    children = get(gca, 'children');
    delete(children(1));
end
close(vO)
close all

