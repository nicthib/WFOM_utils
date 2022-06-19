function [IDX_x, IDX_y, IDX_s] = cnt_silhouette(cnt,IDX,networks,display)
close all
k = size(cnt,1);
n = size(cnt,2);
allc = 1:n;
IDX_x = zeros([size(IDX),k]);
IDX_y = zeros([size(IDX),k]);
IDX_s = zeros([size(IDX),k]);
networks = [networks networks];
cnt_ds = zeros(numel(networks),numel(networks),k);
% Get network means
for ks = 1:k
    for i = 1:numel(networks)
        for j = 1:numel(networks)
            netidx1 = 1:networks(1);
            if i > 1
                netidx1 = [1:networks(i)]+sum(networks(1:i-1));
            end
            netidx2 = 1:networks(1);
            if j > 1
                netidx2 = [1:networks(j)]+sum(networks(1:j-1));
            end
            cnt_ds(i,j,ks) = mean(mean(cnt(ks,netidx1,netidx2)));
        end
    end
end

for i = 1:k
    IDX_xtmp = zeros(size(IDX));
    IDX_ytmp = zeros(size(IDX));
    IDX_stmp = zeros(size(IDX));
    for j = 1:numel(networks)
        netidx = 1:networks(1);
        if j > 1
            netidx = [1:networks(j)]+sum(networks(1:j-1));
        end
            
        nonnetidx = allc(find(~ismember(allc,netidx)));
        x = mean(mean(cnt(i,netidx,netidx)));
        y = mean(mean(cnt(i,nonnetidx,nonnetidx)));
        pixplace = ismember(IDX,netidx);
        IDX_xtmp(pixplace) = x;
        IDX_ytmp(pixplace) = y;
        IDX_stmp(pixplace) = (x-y)/max([x y]);
    end
    IDX_x(:,:,i) = IDX_xtmp;
    IDX_y(:,:,i) = IDX_ytmp;
    IDX_s(:,:,i) = IDX_stmp;
end

% display
if display
IDXn = zeros(size(IDX));
IDXc = zeros(numel(networks),2);
for i = 1:numel(networks)
    netidx = 1:networks(1);
    if i > 1
        netidx = [1:networks(i)]+sum(networks(1:i-1));
    end
    tmp = ismember(IDX,netidx);
    
    IDXn(tmp) = i;
    IDXc(i,:) = regionprops(bwareaopen(tmp,100),'Centroid').Centroid;
end
net_overlay = imerode(logical(imgradient_nic(IDXn)),ones(1));

gfq = @(myStruct,sz) imresize(myStruct.cdata,[sz sz]);

colormap jet
h = tight_subplot(7,3);
for k = 1:7
    axes(h(1+(k-1)*3))
    im = imagesc(IDX_x(:,:,k)); caxis([0 1]); axis image; axis off
    im.AlphaData = ~IDX==0;
    imagesc(imoverlay(gfq(getframe(gca),256),net_overlay,'k'))
    axis image; axis off; caxis([0 1])

    axes(h(2+(k-1)*3))
    im = imagesc(IDX_y(:,:,k)); caxis([0 1]); axis image; axis off
    im.AlphaData = ~IDX==0;
    imagesc(imoverlay(gfq(getframe(gca),256),net_overlay,'k'))
    axis image; axis off; caxis([0 1])

%     axes(h(2+(k-1)*3)); cla
%     im = imagesc(imoverlay(zeros(256,256,3),net_overlay,'w')); axis image; axis off
%     im.AlphaData = ~IDX==0;
% 
%     cmap = jet(100);
%     for i = 1:numel(networks)
%         for j = 1:numel(networks)
%             crow = round(cnt_ds(i,j,k)*100);
%             line([IDXc(i,1) IDXc(j,1)],[IDXc(i,2) IDXc(j,2)],'Color',[cmap(crow,:) .5],'LineWidth',1.5)
%             hold on
%         end
%     end

    axes(h(3+(k-1)*3))
    im = imagesc(IDX_s(:,:,k)); caxis([-1 1]); axis image; axis off
    im.AlphaData = ~IDX==0;
    imagesc(imoverlay(gfq(getframe(gca),256),net_overlay,'k'))
    axis image; axis off; caxis([-1 1]);
end
end
    