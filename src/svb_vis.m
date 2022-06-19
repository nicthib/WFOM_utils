%[rotf,pupil,whisk] = svb_vis(runname,st,c,divs,mapidx)

function [r,p,w] = svb_vis(runname,st,c,divs,mapidx,stimtimes)
close all
load('TF.mat','cm125')
[stb,wz,pz,rz,w,p,r] = getbehaviorz(runname);
w = w-min(w); w = w/max(w);
p = p-min(p); p = p/max(p);
r = r-min(r); r = r/max(r);
k = size(c,1); r = sqrt(size(c,2));

% Visualize
figure('Position',[1000 0 900 1000])
tiledlayout(5,1)
pan xon
lim = [0 11980]/20;
t = (1:11980)/20;

ax(1) = nexttile;
colormap(jet)
imagesc(t,1,(st)'); caxis([1 k])
title('Neural State')
xlim(lim); grid on

% Hidden legend
hold on;
legend_str = {};
for i = 1:k; hidden_h(i) = surf(uint8(i+1-[1 1;1 1]), 'edgecolor', 'none'); end
for i = 1:k; legend_str{end+1} = mat2str(i); end
hold off
uistack(hidden_h, 'bottom');
legend(hidden_h, legend_str,'Location','EastOutside')
yticks([])

ax(2) = nexttile;
imagesc(t,1,(stb{1}+1)); caxis([1 k])
title('Behavioral State')
xlim(lim); grid on
xticks([]); yticks([])

ax(3) = nexttile;
plot(t,w); hold on; plot(t,wz)
title('Whisker motion')
xlim(lim); grid on
ylim([-.1 1.1])
xticks([]); yticks([])

ax(4) = nexttile;
plot(t,r); hold on; plot(t,rz)
title('Locomotion speed')
xlim(lim); grid on
ylim([-.1 1.1])
if ~isempty(stimtimes)
   for i = 1:numel(stimtimes)
       plot([stimtimes(i) stimtimes(i)],[-2 2],'k')
   end
end
xticks([]); yticks([])

ax(5) = nexttile;
plot(t,p); hold on; plot(t,pz)
title('Pupil diameter')
xlim(lim); grid on
ylim([-.1 1.1])
xticks([]); yticks([])

linkaxes(ax,'x')
sgtitle(strrep(runname(1:end-5),'_',' '))

% Behavior summary
cmap = jet(k)*0; cmap2 = jet(k);
figure('Position',[0 0 1500 500])
for j = 1:k
    subplot(2,k,j)
    tmp = divcorrmap(reshape(c(j,:),[r r]),divs);
    im = imagesc(tmp);
    im.AlphaData = ~isnan(tmp);
    caxis([0 1]);
    axis image; axis off
    title(['State ' num2str(j)])
    colormap(jet)
end

for j = 1:k
    subplot(2,k,j+k)
    imagesc(corrIDXmap(cm125.IDXs_b{mapidx},reshape(c(j,:),[r r])))
    caxis([0 1]);
    axis image
    axis off
end
sgtitle(strrep(runname(1:end-5),'_',' '))

% 
% for i = 1:k
%     idx = find(st == i);
%     subplot(6,k,2*k+1:3*k)
%     title('Locomotion')
%     hold on; grid on
%     tmp = r(idx);
%     errorbar(i,mean(tmp),std(tmp),'Color',cmap(i,:))
%     scatter(i,mean(tmp),[],[0 0 0],'filled')
%     set(gca,'xtick',[],'XColor','none')
%     xlim([.6 k+.4])
% 
%     subplot(6,k,3*k+1:4*k)
%     title('Whisking')
%     hold on; grid on
%     tmp = w(idx);
%     errorbar(i,mean(tmp),std(tmp),'Color',cmap(i,:))
%     scatter(i,mean(tmp),[],[0 0 0],'filled')
%     set(gca,'xtick',[],'XColor','none')
%     xlim([.6 k+.4])
%     
%     subplot(6,k,4*k+1:5*k)
%     title('Pupil Diameter')
%     hold on; grid on
%     tmp = p(idx);
%     errorbar(i,mean(tmp),std(tmp),'Color',cmap(i,:))
%     scatter(i,mean(tmp),[],[0 0 0],'filled')
%     set(gca,'xtick',[],'XColor','none')
%     xlabel('State #')
%     xlim([.6 k+.4])
% end
% 
% % Scatter
% for i = 1:k
%     eval(['ax' num2str(i-1) '= subplot(5,k,4*k+i);'])
%     idx = find(st == i);
%     hold on; 
%     tmp1 = p(idx)'; tmp2 = w(idx)';
%     X = [tmp1;tmp2];
%     [x1,x2] = meshgrid(15:40,0:30);
%     x1 = x1(:); x2 = x2(:);
%     xi = [x1(:) x2(:)];
%     ksdensity(X,xi,'plotfcn','contour');
%     axis image; axis off
%     xlabel('Pupil Diameter')
%     ylabel('Whisking')
% end

