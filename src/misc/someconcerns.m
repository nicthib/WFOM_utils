% Do we trust fSVD is even doing the same as PCA?
% Lets test.
% PCA vs SVD
m = makem;
m.loadpct = [0 1];
m.outputs = 'l';
m.PCAcomps = 500;
m.BW = ones(256,256);
m.dsf = 2;
m.nrot = 2;
fulldir = findmousefolder('cm125_3','runB',1);
[~,data1] = LoadData_v21(fulldir,m);

% load raw, then use standard svd
m=makem;
m.loadpct = [0 1];
m.outputs = 'l';
m.BW = ones(256,256);
m.dsf = 2;
m.nrot = 2;
fulldir = findmousefolder('cm125_3','runB',1);
[~,data2] = LoadData_v21(fulldir,m);
[C,S] = pca(reshape(data2.lime,[256^2,size(data2.lime,3)]));
data2.lime = reshape(S(:,1:500)*C(:,1:500)',[256 256 size(data2.lime,3)]);

% measure correlation of entire dataset
similarity = corr(data1.lime(:),data2.lime(:));
pixcorr = zeros(256,256);
for i = 1:256
    for j = 1:256
        pixcorr(i,j) = corr(squeeze(data1.lime(i,j,:)),squeeze(data2.lime(i,j,:)));
    end
end

%% imregister

fixed = data.lime(:,:,1);
[optimizer, metric] = imregconfig('monomodal');
tforms = zeros(3,3,1000);
for i = 1:1000
    TF = imregister_nic(data.lime(:,:,i),fixed,'affine',optimizer,metric);
    tforms(:,:,i) = TF.tform.T;
end

%% To fix the 'DLC' variable issue
vars = who;
for i = 1:numel(vars)
    eval([vars{i}(1:end-3) ' = ' vars{i}])
    eval(['clear ' vars{i}])
end
clear i vars

%% Optimal correction factors
session_names = {'cm125_2','cm126_2','cm127_2','cm128_2'};
run_names = {'runD','runC','runD','runD'};
load('/local_mount/space/dingus/1/RS_analysis/analysis_dec2020/params.mat','IDX');
s = 1; % session ID (1-4)
fulldir = findmousefolder(session_names{s},run_names{s},'1');
m = makem;
m.loadpct = [0 .1];
m.outputs = 'rgl';
m.IDX = IDX.(session_names{s}); m.BW = (m.IDX>0);
m.PCAcomps = 100;
m.dsf = 2; m.nrot = 2;
[m,data] = LoadData_v2(fulldir,m);

% Lowpass
ss = size(data.green);
for i = 1:m.nLEDs
    tmpLED = reshape(data.(m.LEDs{i}),[prod(ss(1:2)) ss(3)]);
    tmpLED(isnan(tmpLED(:,1)),:) = [];
    data.(m.LEDs{i}) = lowpass(tmpLED',5,m.framerate/m.nLEDs);
end

% convert in a loop
Dr = 0:.1:2; Dg = 0:.1:2;
for r = 1:numel(Dr)
    for g = 1:numel(Dg)
        jrgeco = data.lime./((data.red.^Dr(r)).*(data.green.^Dg(g)));
        F0 = mean(jrgeco(:,:,m.baseline),3);
        jrgeco_all{r,g} = jrgeco./repmat(F0,[1 1 size(jrgeco,3)])-1;
    end
end

%
v = VideoWriter(session_names{s});
v.FrameRate = 20;
open(v)
addpath('/local_mount/space/juno/1/Nic/workflow_nic/utils/')
r = [1 1 20 8]; g = [1 20 1 4]; % Dr and Dg values to show x10
close all; figure('Position',[0 0 1000 500])
ha = tight_subplot(2,numel(r),[.01 .01],[0 .1],[.01 .01]);
sgtitle(strrep(session_names{s},'_',' '))
for i = 1:100
    for j = 1:numel(r)
        axes(ha(j))
        frame = jrgeco_all{r(j),g(j)}(:,:,i);
        imagesc(frame)
        title(['Dr = ' mat2str(Dr(r(j))) ', Dg = ' mat2str(Dg(g(j)))])
        caxis([-.1 .1]);
        colormap gray; axis image; axis off;
        axes(ha(j+numel(r)))
        imagesc(imgradient(frame))
        caxis([0 .1])
        colormap gray; axis image; axis off;
    end
    f = getframe(gcf);
    writeVideo(v,f);
end
close(v)
close all

%% Filter comparison (1/15/21)
savedir = '/local_mount/space/chronos/2/cmdata_analysis/RS_H_1_15_21';
DLCdir = '/local_mount/space/dingus/1/RS_analysis/DLC/final';
load(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_jan2021','params.mat'))
filter_file = '/local_mount/space/dingus/1/RS_analysis/analysis_jan2021/lowpassFilter_maximallyFlat.mat';
i = 17; j = 1; % cm128_2 run B

file_id = [sessions_for_analysis{i} '_' runs_for_analysis{j}];
load_dir = findmousefolder(sessions_for_analysis{i},runs_for_analysis{j},'1');
m = makem;
m.Dr = 1.1; m.Dg = .2; % Based on recent analysis
m.outputs = 'n'; m.dsf = 2; m.nrot = 2; m.loadpct = [0 .2];
m.IDX = IDX.(sessions_for_analysis{i}); m.BW = (m.IDX>0);
m.bkgsub = bkg;
% First load with no filter
[~,data1] = LoadData_v21(load_dir,m);
% Add filter file and load again
m.filter = filter_file;
[~,data2] = LoadData_v21(load_dir,m);

data3 = cat(2,data1.jrgeco,data2.jrgeco,data2.jrgeco-data1.jrgeco);
%quickvideo(data3,20,[-.1 .1],'Filter_comparison.avi')

%% LED comparison
LED_resid = zeros(256,256*3,1000);
LEDs = {'lime','green','red'};
for i = 1:numel(LEDs)
    LED_resid(:,[1:256]+256*(i-1),:) = data2.(LEDs{i})(:,:,1:1000)-data1.(LEDs{i})(:,:,1:1000);
end
quickvideo(LED_resid,20,[-400 400],'Filter_comparison.avi')

%% TC example
close all
for i = 1:numel(LEDs)
    subplot(3,2,2*(i-1)+1)
    plot(linspace(0,20,401),data1.H.(LEDs{i})(20,100:500)-data1.H.(LEDs{i})(20,100))
    hold on
    plot(linspace(0,20,401),data2.H.(LEDs{i})(20,100:500)-data2.H.(LEDs{i})(20,100))
    plot(linspace(0,20,401),data2.H.(LEDs{i})(20,100:500)-data1.H.(LEDs{i})(20,100:500))
    legend('Raw','Filtered','Residual')
    xlabel('Time (sec)')
    title(LEDs{i})
    ylim([-100 100])
    
    subplot(3,2,i*2)
    
    [~,f1,Spectrum1,~] = fastFourTrans(data1.H.(LEDs{i})(20,:),20,20);
    [~,f2,Spectrum2,~] = fastFourTrans(data2.H.(LEDs{i})(20,:),20,20);
    plot(f1,Spectrum1); hold on
    plot(f2,Spectrum2)
    ylim([0 10])
    legend('raw','filtered')
    xlabel('Freq (hz)')
end

%% TC example
close all
plot(linspace(0,20,401),data1.H.jrgeco(20,100:500))
hold on
plot(linspace(0,20,401),data2.H.jrgeco(20,100:500))
plot(linspace(0,20,401),H1(20,100:500))
legend('Raw','Filtered','Residual')
xlabel('Time (sec)')
ylabel('dF/F')

%% Stationary FC
subplot(131)
imagesc(corr(data1.H.jrgeco'))
axis image; axis off; caxis([0 1]); title('Raw')
colorbar
subplot(132)
imagesc(corr(data2.H.jrgeco'))
axis image; axis off; caxis([0 1]); title('Filter')
colorbar
subplot(133)
imagesc(corr(data2.H.jrgeco')-corr(data1.H.jrgeco'))
axis image; axis off; title('Difference')
caxis([-.01 .01])
colorbar
colormap jet

%% Loading raw spectra files and processing them
clear
filename = 'spectra_010621_WFOM2';
path_name = '/local_mount/space/dingus/1/RS_analysis/spectra/spectra_010621_WFOM2';
files = dir(path_name);
files = files(3:end);
for i = 1:numel(files)
    id{i} = files(i).name(strfind(files(i).name,'_')+1:strfind(files(i).name,'.ProcSpec')-1);
    eval(['spectra_' id{i} '=importdata(''' fullfile(files(i).folder,files(i).name) ''');']);
end

% temporary fix for missing blue spectra
LEDs = {'lime','green','red','blue'};
for i = 1:numel(LEDs)
    var = ['spectra_' LEDs{i}];
    eval([var '(:,2) =' var '(:,2)-spectra_dark(:,2);']);
end

plotcolors = {[0 0 1],[0 0 0],[.5  0 0],[0 1 0],[0 .5 0],[1 0 0]};

close all
for i = 1:numel(id)
   var = ['spectra_' id{i}];
   eval(['plot(' var '(:,1),' var '(:,2),''Color'',' mat2str(plotcolors{i}) ')']) 
   hold on
end
legend(id)
ylim([0 60000])
xlim([500 700])

save([filename '.mat'],'spectra_red','spectra_green','spectra_blue')

%% Comparing old spectra to new for heomdynamics

clear
load(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_jan2021','params.mat'))
filter_file = '/local_mount/space/dingus/1/RS_analysis/analysis_jan2021/lowpassFilter_maximallyFlat.mat';
spectra_file_old = '/local_mount/space/juno/1/Software/MIAO/Spectra/LED_spectra_0711.mat';
spectra_file_new = '/local_mount/space/dingus/1/RS_analysis/spectra/spectra_010621_WFOM2.mat';
i = 3; j = 7;
load_dir = findmousefolder(sessions_for_analysis{i},runs_for_analysis{j},'1');

m = makem;
m.outputs = 'odn'; m.dsf = 2; m.nrot = 2; m.loadpct = [0 .1];%.9995];
m.filter_file = filter_file;
m.IDX = IDX.(sessions_for_analysis{i}); m.BW = (m.IDX > 0);
m.bkgsub = bkg;
m.spectra_file = spectra_file_old;
m.baseline = 100:200;
[m1,data1] = LoadData_v21(load_dir,m);
m.spectra_file = spectra_file_new;
[m2,data2] = LoadData_v21(load_dir,m);
epoch = 200:1000;

%%
chbo1 = data1.chbo(:,:,epoch);
chbo2 = data2.chbo(:,:,epoch);
chbr1 = data1.chbr(:,:,epoch);
chbr2 = data2.chbr(:,:,epoch);
chbt1 = chbo1+chbr1;
chbt2 = chbo2+chbr2;

jrgeco = data1.jrgeco(:,:,epoch);

%%
close all
figure('Position',[0 0 1000 1000])
subplot(3,3,8:9)
t = linspace(0,numel(epoch)/20,numel(epoch));
yyaxis left
hold on
plot(t,data1.H.chbo(31,epoch),'r-')
plot(t,data1.H.chbr(31,epoch),'b-')
plot(t,data1.H.chbo(31,epoch)+data1.H.chbr(31,epoch),'g-')

ylim([-10 10]*3e-6)
ylabel('Hb')
yyaxis right
plot(t,data1.H.jrgeco(31,epoch),'k')
ylim([-.2 .2])
ylabel('jrgeco')
xlabel('time(sec)')
l = line([0 0],[-.2 .2],'Color','k');
v = VideoWriter('Hb_comparison.avi');
v.FrameRate = 10;
open(v)
for i = 1:numel(epoch)
    subplot(331)
    imagesc(chbo1(:,:,i))
    axis image; axis off; caxis([-10 10]*2e-6)
    title('CHBO (old)')
    subplot(332)
    imagesc(chbr1(:,:,i))
    axis image; axis off; caxis([-10 10]*2e-6)
    title('CHBR (old)')
    subplot(333)
    imagesc(chbt1(:,:,i))
    axis image; axis off; caxis([-10 10]*2e-6)
    title('CHBT (old)')
    subplot(334)
    imagesc(chbo2(:,:,i))
    axis image; axis off; caxis([-10 10]*2e-6)
    title('CHBO (new)')
    subplot(335)
    imagesc(chbr2(:,:,i))
    axis image; axis off; caxis([-10 10]*2e-6)
    title('CHBR (new)')
    subplot(336)
    imagesc(chbt2(:,:,i))
    axis image; axis off; caxis([-10 10]*2e-6)
    title('CHBT (new)')
    subplot(337)
    imagesc(jrgeco(:,:,i))
    axis image; axis off;
    title('jRGECO')
    colormap jet
    caxis([-.1 .1])
    
    subplot(3,3,8:9)
    xlim([-5 5]+i/20)
    l.XData = [i i]/20;
    f = getframe(gcf);
    writeVideo(v,f);
end
close(v)

%%
close all
n = [10 24 31 2 20];
t = linspace(0,numel(epoch)/20,numel(epoch));
figure
for i = 1:numel(n)
    subplot(5,3,[1:2]+(i-1)*3)
    yyaxis left
    plot(t,data1.H.jrgeco(n(i),epoch),'k')
    xlim([0 40])
    yyaxis right
    hold on
    plot(t,(data1.H.chbo(n(i),epoch)),'r-')
    plot(t,(data2.H.chbo(n(i),epoch)),'r--')
    plot(t,(data1.H.chbr(n(i),epoch)),'b-')
    plot(t,(data2.H.chbr(n(i),epoch)),'b--')
    plot(t,(data1.H.chbo(n(i),epoch)+data1.H.chbr(n(i),epoch)),'g-')
    plot(t,(data2.H.chbo(n(i),epoch)+data2.H.chbr(n(i),epoch)),'g--')
    
    subplot(5,3,3+(i-1)*3)
    imshowpair(logical(imgradient_nic(m.IDX)),m.IDX==n(i)); axis image; axis off
end

%% Filter LEDs vs. Filter converted data
clear
load(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_jan2021','params.mat'))
filter_file = '/local_mount/space/dingus/1/RS_analysis/analysis_jan2021/lowpassFilter_maximallyFlat.mat';
spectra_file = '/local_mount/space/dingus/1/RS_analysis/spectra/spectra_010621_WFOM2.mat';
i = 3; j = 7;
load_dir = findmousefolder(sessions_for_analysis{i},runs_for_analysis{j},'1');

m = makem;
m.outputs = 'odn'; m.dsf = 2; m.nrot = 2; m.loadpct = [0 .1];%.9995];
m.IDX = IDX.(sessions_for_analysis{i}); m.BW = (m.IDX > 0);
m.bkgsub = bkg;
m.spectra_file = spectra_file;
m.baseline = 100:200;
[m1,data1] = LoadData_v21(load_dir,m);
m.filter_file = filter_file;
[m2,data2] = LoadData_v21(load_dir,m);
epoch = 200:1000;

%% 
load(m.filter_file,'Hd')
delay = mean(grpdelay(Hd));
tmp = filter(Hd,data1.H.chbo(20,epoch),2);

figure
plot(data1.H.chbo(20,epoch))
hold on
plot(data2.H.chbo(20,epoch))
plot(tmp(delay:end))
legend('Unfiltered','Filter LEDs','Filter converted')

%% Comparing denosied vs noised data
load(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_dec2020','params.mat'))
m = makem;
m.outputs = 'odn';
m.dsf = 2;
m.nrot = 2;
m.loadpct = [0 1];
m.BW = ones(256,256);
mouse = sessions_for_analysis{1};
m.IDX = IDX.(mouse);
fulldir = findmousefolder(mouse,'runB','1');
[~,data1] = LoadData_v21(fulldir,m);
m.PCAcomps = 500;
[m,data2] = LoadData_v21(fulldir,m);

%% Calculate silhouette
silh1 = IDX_H_std(m.IDX,reshape(data1.jrgeco(:,:,5400:6200),[256^2 801]));
silh2 = IDX_H_std(m.IDX,reshape(data2.jrgeco(:,:,5400:6200),[256^2 801]));

%% Show silh
subplot(121)
imagesc(silh1)
caxis([0 1])
axis image
axis off
title('Silhouette - noisy data')
colorbar
subplot(122)
imagesc(silh2)
caxis([0 1])
axis image
axis off
colormap jet
title('Silhouette - denoised data')
colorbar

%% Show corrs
HcorrIDX = zeros(256,256);
for i = 1:92
    Hcorr(i) = corr(data1.H.jrgeco(i,1:1000)',data2.H.jrgeco(i,1:1000)');
    HcorrIDX(m.IDX==i) = Hcorr(i);
end




