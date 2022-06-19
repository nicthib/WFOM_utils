%% Possible new autocrop?
addpath('/local_mount/space/chronos/2/RS_SVD_2020')
%load('New_BWs.mat'); BWs = BW; clear BW
clearvars -EXCEPT IDX_final H_final runnames exemplars BWs
runname = 'cm125_2_runD.mat';
load(runname)
rotf = getrotf(m);
m.baseline = baselinefromrot(rotf,100,1);
for d = 1:m.nLEDs
    data.(m.LEDs{d}) = convertSVD(C.(m.LEDs{d}),S.(m.LEDs{d}),1:500,m.nanidx,256);
end
data.jrgeco = data.lime./((data.red.^m.Dr).*(data.green.^m.Dg));
bgGG = mean(data.jrgeco(:,:,m.baseline),3);
data.jrgeco = data.jrgeco./repmat(bgGG,[1 1 size(data.jrgeco,3)])-1;
data.jrgeco = reshape(data.jrgeco,[256^2,size(data.jrgeco,3)]);
jrgeco = data.jrgeco(:,1:11000);
jrgeco(isinf(jrgeco)) = 0; jrgeco(isnan(jrgeco)) = 0;
rotf = rotf(1:11000);
clear data

% hp filter
jrgeco_filt = zeros(size(jrgeco))';
fs = m.framerate/3;
f_n = [0.5, 0.6]*2/fs; % stopband in normalized frame rate, so here it is 0.5Hz hp-filter
a = [1 0];
b_hp = firpm(500,[0 f_n 1], [0 0 1 1]); % filter order is 500
non_nan_ind = ~isnan(jrgeco(:,1));
jrgeco_filt(:,non_nan_ind) = filtfilt(b_hp,a,jrgeco(non_nan_ind,:)');
jrgeco_filt = jrgeco_filt';
jrgeco = jrgeco_filt; clearvars -EXCEPT jrgeco m rotf IDX_final

%%
kmeans_idx = baselinefromrot(rotf,500,0); kmeans_idx = kmeans_idx(50:450);
%moving_idx = baselinefromrot(-smooth(rotf,50),300,-.2);
%exemplar = jrgeco_filt(:,kmeans_idx);
%exemplar_moving = jrgeco_filt(:,moving_idx);
% for i = 1:3
%     IDX = kmeans(exemplar,100,'Distance','Correlation');
%     IDX = reshape(IDX,[256 256]);
%     BW = imfill(logical(cleanupIDX(IDX,200)),'holes');
%     exemplar(BW(:)==0,:) = NaN;
% end
% BW = imerode(bwconvhull(bwareaopen(BW,10000)),ones(5,5));
BW = BWs.(m.mouse);
exemplar = jrgeco_filt(:,kmeans_idx);
exemplar(BW(:)==0,:) = NaN;

%
imagesc(reshape(std(exemplar(:,1:100),[],2),[256 256]));
[BW_L,xi,yi] = roipoly;
mid_vec = [-xi(1)+xi(2),yi(1)-yi(2),0];
BW_L = and(BW,BW_L);
BW_R = and(BW,~BW_L);
close all
exemplar_L = exemplar;
exemplar_R = exemplar;
exemplar_L(BW_L(:)==0) = NaN;
exemplar_R(BW_R(:)==0) = NaN;

IDX_L = kmeans(exemplar_L,50,'Distance','Correlation','Replicates',10);
IDX_L = reshape(IDX_L,[256 256]);
H = getHfromKmeans(exemplar_L,IDX_L,0);
IDX_R = kmeans(exemplar_R,50,'Distance','Correlation','Start',H);
IDX_R = reshape(IDX_R,[256 256]);
IDX_LR = sum(cat(3,IDX_L,IDX_R),3,'omitnan');
H_final.(m.mouse) = getHfromKmeans(exemplar,IDX_LR,0);
IDX_final.(m.mouse).LR = IDX_LR;
IDX_final.(m.mouse).L = IDX_L;
IDX_final.(m.mouse).R = IDX_R;
IDX_final.(m.mouse).BW_L = BW_L;
IDX_final.(m.mouse).BW_R = BW_R;
runnames.(m.mouse) = runname;
exemplars.(m.mouse).RS = exemplar;
%exemplars.(m.mouse).AS = exemplar_moving;

clearvars -]EXCEPT IDX_final H_final runnames exemplars BWs

%% Another approach
close all
figure('Position',[0 0 1000 500])
IDX_LR_refined = zeros(256,256);
currcomp = zeros(256,256);
i = min(find(~ismember(1:50,unique(IDX_LR_refined(:)))));
subplot(131); imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
title('Press Up arrow to store component.')
subplot(133); imagesc(IDX_LR_refined); axis image; axis off
while true
    [x,y] = ginput(1); x = round(x); y = round(y);
    if BW_L(y,x) == 1
        currcomp(IDX_L==IDX_L(y,x)) = ~currcomp(y,x)*2;
    else
        currcomp(IDX_R==IDX_R(y,x)) = ~currcomp(y,x)*2;
    end
    subplot(131); imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
    
    subplot(132);
    if sum(currcomp(:)) > 0
        c = corr(H',getHfromKmeans(exemplar,currcomp/2,0)');
        cIDX = zeros(256,256);
        for j = 1:max(IDX_LR(:))
            cIDX(IDX_LR==j) = c(j);
        end
        imagesc(cIDX); axis image; axis off; caxis([0 1])
    end
    title('Press Up arrow to store component.')
    k = waitforbuttonpress;
    v = double(get(gcf,'CurrentCharacter'));
    if v == 30 % U arrow (STORE)
        IDX_LR_refined = IDX_LR_refined + (currcomp/2).*i;
        subplot(133); imagesc(IDX_LR_refined); axis image; axis off
        currcomp = zeros(256,256);
        i = min(find(~ismember(1:50,unique(IDX_LR_refined(:)))));
    end
end

%%
IDX_LR_refined = sortkmeans(IDX_LR_refined);
IDX_final.(m.mouse) = roundIDX(IDX_LR_refined .* bwareaopen(logical(IDX_LR_refined),1000));

%%
IDXs = fieldnames(IDX_final);
IDXbest = IDX_final.cm127_2;
for i = 1:max(IDXbest(:))
    subplot(144)
    imshowpair(logical(imgradient(IDXbest)),IDXbest==i); axis image; axis off
    IDXtmp = IDXbest==i;
    for j = 1:3
        IDXtest = IDX_final.(IDXs{j});
        overlay = IDXtest.*IDXtmp;
        overlay(overlay==0) = NaN;
        bestcomp = mode(overlay(:));
        subplot(1,4,j)
        imshowpair(logical(imgradient(IDXtest)),IDXtest==bestcomp); axis image; axis off
    end
    pause
end

%% Approach 1 but automatic
criteria = {'Area','MajorAxisLength','MinorAxisLength','Centroid'};
for i = 1:max(IDX_L(:))
    tmpL = bwareaopen(IDX_L==i,100);
    tmpR = bwareaopen(IDX_R==i,100);
    Lp = regionprops(tmpL,criteria);
    Rp = regionprops(tmpR,criteria);
    for j = 1:numel(criteria)
        try
            sim_valL = 0; sim_valR = 0;
            for k = 1:numel(Lp)
                sim_valL = sim_valL + Lp(k).(criteria{j});
            end
            for k = 1:numel(Rp)
                sim_valR = sim_valR + Rp(k).(criteria{j});
            end
            if j ~= 1
                sim_valL = sim_valL/numel(Lp);
                sim_valR = sim_valR/numel(Rp);
            end
            if j == 4
                vec_L = [128-sim_valL(1),sim_valL(2)-128,0]; % Rise/run from center
                vec_R = [128-sim_valR(1),sim_valR(2)-128,0]; % Rise/-run from center
                x_L = cross(mid_vec,vec_L); x_R = cross(mid_vec,vec_R);
                sim_val(i,j) = abs(x_L(3)/x_R(3));
            else
                sim_val(i,j) = sim_valL/sim_valR;
            end
            if sim_val(i,j) > 1
                sim_val(i,j) = sim_val(i,j).^-1;
            end
        catch
        end
    end
end

quality_tmp = zeros(256,256);
for i = 1:max(IDX_LR(:))
   quality_tmp(IDX_LR==i)=mean(sim_val(i,:));
end

good_comps = find(mean(sim_val,2)>=.8);
bad_comps = 1:50; bad_comps(good_comps) = [];
IDX_LR_refined = IDX_LR;
IDX_LR_refined(ismember(IDX_LR,bad_comps)) = 0;

