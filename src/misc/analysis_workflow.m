%% Initialize
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
clearvars -EXCEPT pths
%% Load
m = makem; try m.BW = BW; catch; end
m.mouse = 'cm128_4';  m.dsf = 2; m.nrot = 2; %m.PCAcomps = 250; 
m.run = 'runC'; m.stim = '1'; m.ccddir = findmousefolder(m.mouse);
m.smooth = []; m.loadpct = [0 .3]; %m.baseline = 6000:7000;
fulldir = findmousefolder(m.mouse,m.run,m.stim); m.outputs = 'rglodn';
[m,data] = LoadData_v2(fulldir,m);
jrgeco = data.jrgeco; chbo = data.chbo; chbr = data.chbr; %clear data
ss = size(jrgeco); m.sz = ss(1);

%% Obtaining H_n and W_n
ncomps = 100; BW = m.BW;
BW = double(BW); BW(BW==0) = NaN;
[IDX,C] = kmeans(reshape(jrgeco(:,:,1000:1999),[prod(ss(1:2)) 1000]),ncomps,'MaxIter',100,'Display','off','OnlinePhase','off','Distance','correlation');
close all; imagesc(reshape(IDX,[m.sz m.sz]));

% Get basis TC's, perform NNLSQ
H_n = getHfromKmeans(jrgeco,IDX,m,0);
H_n(isnan(H_n(:,1)),:) = []; % remove empty components
[~,W_n] = LSQanalysis(jrgeco(:,:,m.baseline),H_n(:,m.baseline));

% Arrange comps
I = arrangecomps(W_n); % arranges kmeans components from top to bottom
H_n = H_n(I,:); W_n = W_n(:,I);
cmap = jet(size(H_n,1));

%% Make videos of model vs raw
clear im
cax = [0 .1];
savedir = pths.video;
f_num = ss(3);
for v = 1 % makes ~90 second chunk videos
    filename = [m.mouse '_' m.run '_' mat2str(v) '_neural.avi'];
    vidObj = VideoWriter(fullfile(savedir,filename));
    vidObj.FrameRate = m.framerate/3;
    e = [1:f_num] + (v-1)*f_num;
    open(vidObj)
    parfor i = 1:numel(e)
        d1 = padarray(im2u8sc(reshape(W_n*diag(H_n(:,e(i))+.02)*cmap,[m.sz m.sz 3]),cax),[20 20],'both');
        d0 = padarray(im2u8sc(repmat(jrgeco(:,:,e(i))+.04,[1 1 3]),cax),[20 20],'both');
        d0 = insertText(d0,[10 10],sprintf('%0.1f sec',round(e(i)*100/(m.framerate/3))/100),'FontSize',15,'TextColor','w','BoxOpacity',0);
        tmp = permute(reshape(flipud(H_n(:,1:i)),[size(H_n,1) 1 i]).*repmat(cmap,[1 1 i]),[1 3 2]);
        im(i).cdata = [imresize([d0 d1],[360 720]*2); im2u8sc(imresize(tmp,[400 1440]),cax)];
        im(i).colormap = [];
    end
    disp(['Writing file ' filename])
    writeVideo(vidObj,im);
    close(vidObj)
    clear im
end
disp('Done')

%% Make FLIR video
f_start = (m.spoolsLoaded(1)-1)*m.numFramesPerSpool;
f_num = (m.spoolsLoaded(end)*m.numFramesPerSpool-(m.spoolsLoaded(1)-1)*m.numFramesPerSpool)/m.nLEDs;
BatchConvertFLIR(m.mouse,m.run,f_start,f_num,pths.video)

%% Show NNLSQ output
close all
figure('Position',[0,0,700,700]);
for i = 1:size(W_n,2)
    subplot(8,8,i)
    imagesc(reshape(W_n(:,i)*cmap(i,:),[m.sz m.sz 3]))
    title(mat2str(i))
    axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end

%% Further analysis
%% generate correlation matrix
H_corr = temporal_corr(H_n,60);
H_corr = permute(H_corr,[3 1 2]);

%% kmeans correlations
[IDX2, C2, SUMD] = kmeans(reshape(H_corr,[size(H_corr,1),size(H_corr,2)^2]), 2);
figure; plot(IDX2)
figure; imagesc(IDX2')

%% 
close all
subplot(211)
imagesc(IDX2')
subplot(212)
plot(diff(smooth(double(unwrap(m.aux(2,:),100)))))
xlim([0 numel(m.aux(2,:))])

%% Silhouette
ncomps = [10 20 30 40 50 60 70 80 90];
for i = 6
    nc = ncomps(i); BW = m.BW;
    BW = double(BW); BW(BW==0) = NaN;
    [IDX] = kmeans(reshape(jrgeco(:,:,2000:2999),[prod(ss(1:2)) 1000]),nc,'MaxIter',100,'Display','off','OnlinePhase','off','Distance','correlation');
    subplot(3,3,i)
    [S,H] = silhouette(reshape(jrgeco(:,:,1:1000),[prod(ss(1:2)) 1000]),IDX,'correlation');
    title(nanmean(S))
end

