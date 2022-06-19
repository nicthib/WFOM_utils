clear
filename = 'cm126_3_runB.mat';
load(filename);
rotf = getrotf(m);
m.baseline = baselinefromrot(rotf,100,1);
disp(['Converting ' filename])
for d = 1:m.nLEDs
    data.(m.LEDs{d}) = convertSVD(C.(m.LEDs{d}),S.(m.LEDs{d}),1:500,m.nanidx,256);
end
data.jrgeco = data.lime./((data.red.^m.Dr).*(data.green.^m.Dg));
[data.chbo,data.chbr,~] = convert_mariel_MIAO(data.green,data.red,'g','r',m.baseline,m.greenfilter);
bgGG = mean(data.jrgeco(:,:,m.baseline),3);
data.jrgeco = data.jrgeco./repmat(bgGG,[1 1 size(data.jrgeco,3)])-1;
% reshape to speed up H


%%
close all
samp = squeeze(mean(data.jrgeco(:,:,100:200),3));
imagesc(samp); axis image; colormap gray;
title('Click bottom then top')
[x,y] = ginput(2); x=round(x); y=round(y);
close all

%% Ring ROI
seed_pts = [];
b_r = sqrt(diff(x)^2+diff(y)^2)/2-10;
% Ring 1
ringn = [12 20 30 30]; 
d = [b_r/3 2*b_r/3 b_r];
for j = 1:numel(d)
    angi = 360/ringn(j);
    for i = 0:ringn(j)/2-1
        seed_pts(end+1,:) = radialpt(x,y,-i*angi,d(j));
        seed_pts(end+1,:) = radialpt(x,y,(i+1)*angi,d(j));
    end
end

% Plot seeds
close all
imagesc(samp)
axis image
axis off
hold on
line(x,y,'Color','w','LineWidth',2)
scatter(seed_pts(:,1),seed_pts(:,2))

% Get H
C = []; sz = 256; kmeans_idx = 1:size(data.jrgeco,3); bx = 10;
data.jrgeco = reshape(data.jrgeco,[sz sz size(data.jrgeco,ndims(data.jrgeco))]);
for i = 1:size(seed_pts,1)
   r = round(seed_pts(i,2)); c = round(seed_pts(i,1));
   C(i,:) = squeeze(nanmean(nanmean(data.jrgeco(r-bx:r+bx,c-bx:c+bx,kmeans_idx))));
end

%% Perform kmeans
IDX_clk = kmeans(reshape(data.jrgeco(:,:,kmeans_idx),[sz^2,numel(kmeans_idx)]),size(seed_pts,1),'Distance','Correlation','Start',C);
IDX_clk = reshape(IDX_clk,[256 256]);

%%
close all
subplot(121)
imagesc(samp)
axis image
axis off
hold on
line(x,y,'Color','w','LineWidth',2)
for i = 1:size(seed_pts,1)
    rectangle('Position',[seed_pts(i,1)-5 seed_pts(i,2)-5 10 10],'Facecolor','w','EdgeColor','None')
    text(seed_pts(i,1),seed_pts(i,2),mat2str(i),'HorizontalAlignment','Center')
end

subplot(122)
imagesc(reshape(IDX,[sz sz])); axis image; axis off; colormap jet
hold on
caxis([1 size(seed_pts,1)])
cmap = hsv(size(seed_pts,1));
for i = 1:size(seed_pts)
    tmp = zeros(size(IDX));
    tmp(IDX==i) = 1; tmp = bwareaopen(tmp,100);
    stats = regionprops(tmp,'centroid');
    c = stats.Centroid;
    line([seed_pts(i,1) c(1)],[seed_pts(i,2) c(2)])
    scatter(c(1),c(2),150,'.','k')
    %text(seed_pts(i,1),seed_pts(i,2),mat2str(i),'Color',cmap(i,:),'HorizontalAlignment','Center')
    hold on
end
scatter(seed_pts(:,1),seed_pts(:,2),150,'.','w')

%% plot
close all
subplot(221)
imagesc(samp); axis image; axis off; colormap jet
caxis([0 10000])
hold on
line(x,y,'Color','w','LineWidth',2)
subplot(222)
imagesc(ta)
axis image; axis off; colorbar
subplot(223)
imagesc(ta2)
axis image; axis off; colorbar
subplot(224)
imagesc(ta3)
axis image; axis off; colorbar

figure
subplot(121)
imagesc(ta3)
axis image; axis off; colorbar; colormap jet
subplot(122)
imagesc(reshape(IDX,[sz sz]))
axis image; axis off; colorbar

%%
for i = 1:size(seed_pts)
    diff_seed = IDX(round(seed_pts(i,2))-5:round(seed_pts(i,2))+5,round(seed_pts(i,1))-5:round(seed_pts(i,1))+5)==i;
    IDX_acc(i) = sum(diff_seed(:))/numel(diff_seed(:));
end

%% OLD CODE
%%
Cx = mean(x); Cy = mean(y);
ang = -atan((y(2)-y(1))/(x(2)-x(1)))*180/pi;
ta = zeros(size(samp));
for i = 1:size(samp,1)
    for j = 1:size(samp,2)
        opp = Cy-i; adj = -(Cx-j);
        ta(i,j) = atand(opp/adj);
        if adj < 0
           ta(i,j) = ta(i,j)+180;
        elseif opp < 0
            ta(i,j) = ta(i,j)+360;
        end
    end
end
ta=ta-ang;
for i = 1:numel(ta(:))
    if ta(i) < 0
        ta(i) = ta(i)+360;
    elseif ta(i) > 360
        ta(i) = ta(i) - 360;
    end
end
npie = 20;
ta = ta-mod(ta,360/npie);
ta=ta/(360/npie)+1;
%%
ta2 = zeros(size(samp));
for i = 1:npie
    if i <= npie/2
        ta2(ta==i) = i*2;
    else
        ta2(ta==i) = (-i+npie)*2+1;
    end
end
Ds = 4; rad=150;
for i = 1:size(samp,1)
    for j = 1:size(samp,2)
        D = sqrt((Cx-j)^2+(Cy-i)^2);
        if Ds > 1
            ring(i,j) = D-mod(D,rad/Ds);
        end
    end
end
ring = ring/(rad/Ds);
ta2 = (ta2+ring/Ds)*Ds;
ta3 = ta2;
ta3(ta3==0) = NaN;
ta3(samp==0) = NaN;
ta3=ta3-Ds+1;