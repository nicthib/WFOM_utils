function RS_summary_fig(H,m)
% Generate summary figure
% A few TCs

epoch = 100:1000;
figure('Position',[0 0 1920 1080]);
n = round(linspace(1,size(H.jrgeco,1),4));
t = linspace(0,numel(epoch)/20,numel(epoch));
sp_sz = [9 8];
for i = 1:numel(n)
    subplot(sp_sz(1),sp_sz(2),sp_extent([i 1],[1 3],sp_sz))
    yyaxis left
    plot(t,H.jrgeco(n(i),epoch),'k')
    ylabel('jRGECO')
    %xlim([0 40])
    yyaxis right
    ylabel('CHB')
    hold on
    plot(t,(H.chbo(n(i),epoch)),'r-')
    plot(t,(H.chbr(n(i),epoch)),'b-')
    plot(t,(H.chbo(n(i),epoch)+H.chbr(n(i),epoch)),'g-')
    subplot(sp_sz(1),sp_sz(2),sp_extent([i 4],[1 1],sp_sz))
    imshowpair(logical(imgradient_nic(m.IDX)),m.IDX==n(i)); axis image; axis off
end

% Rotary to check w/ TCs
subplot(sp_sz(1),sp_sz(2),sp_extent([5 1],[1 3],sp_sz))
plot(linspace(0,numel(epoch)/(m.framerate/m.nLEDs),numel(epoch)),m.rotf(epoch))
legend('Movement')
xlabel('Time (sec)')

% Baseline images
subplot(sp_sz(1),sp_sz(2),sp_extent([1 5],[2 2],sp_sz))
imshowpair(m.bl_im.lime,logical(imgradient_nic(m.IDX))); axis image; axis off;
title('Lime baseline w/ IDX')
subplot(sp_sz(1),sp_sz(2),sp_extent([1 7],[2 2],sp_sz))
imagesc(m.bl_im.green); axis image; axis off;
title('Green baseline')
subplot(sp_sz(1),sp_sz(2),sp_extent([3 5],[2 2],sp_sz))
imagesc(m.bl_im.red); axis image; axis off;
title('Red baseline')
subplot(sp_sz(1),sp_sz(2),sp_extent([3 7],[2 2],sp_sz))
imagesc(m.jrgeco_quality_im); axis image; axis off;
title('jRGECO quality')
caxis([0 .3]); colorbar

% Spectra
subplot(sp_sz(1),sp_sz(2),sp_extent([7 1],[1 2],sp_sz))
hold on
plot(m.spectra_f,log(m.spectra.lime_prefilt),'k')
plot(m.spectra_f,log(m.spectra.lime_postfilt),'r')
title('Spectra (lime)')
legend('Pre','Post')

subplot(sp_sz(1),sp_sz(2),sp_extent([8 1],[1 2],sp_sz))
hold on
plot(m.spectra_f,log(m.spectra.green_prefilt)+.5,'k')
plot(m.spectra_f,log(m.spectra.green_postfilt)+.5,'r')
title('Spectra (green)')
legend('Pre','Post')

subplot(sp_sz(1),sp_sz(2),sp_extent([9 1],[1 2],sp_sz))
hold on
plot(m.spectra_f,log(m.spectra.red_prefilt)+1,'k')
plot(m.spectra_f,log(m.spectra.red_postfilt)+1,'r')
title('Spectra (red)')
legend('Pre','Post')

% TCs
subplot(sp_sz(1),sp_sz(2),sp_extent([7 3],[1 2],sp_sz))
hold on
plot(t,m.TCs.lime_prefilt(epoch),'k')
plot(t,m.TCs.lime_postfilt(epoch),'r')
title('LED TC (lime)')
legend('Pre','Post')

subplot(sp_sz(1),sp_sz(2),sp_extent([8 3],[1 2],sp_sz))
hold on
plot(t,m.TCs.green_prefilt(epoch),'k')
plot(t,m.TCs.green_postfilt(epoch),'r')
title('LED TC (green)')
legend('Pre','Post')

subplot(sp_sz(1),sp_sz(2),sp_extent([9 3],[1 2],sp_sz))
hold on
plot(t,m.TCs.red_prefilt(epoch),'k')
plot(t,m.TCs.red_postfilt(epoch),'r')
title('LED TC (red)')
legend('Pre','Post')

% Movement
subplot(sp_sz(1),sp_sz(2),sp_extent([6 6],[2 3],sp_sz))
hold on
plot(linspace(0,m.movielength,numel(m.whisk)),m.whisk,'r-')
plot(linspace(0,m.movielength,numel(m.rotf)),zscore(m.rotf),'b-')
plot(linspace(0,m.movielength,numel(m.pupil)),m.pupil,'g-')
xlim([0 m.movielength])
legend('Whisking','Locomotion','Pupil')

subplot(sp_sz(1),sp_sz(2),sp_extent([8 6],[2 3],sp_sz))
t = linspace(0,numel(m.jrgeco_quality_TC)/(m.framerate/m.nLEDs),numel(m.jrgeco_quality_TC));
plot(t,m.jrgeco_quality_TC,'k-')
xlabel('Time (sec)')
xlim([0 m.movielength])
legend('Conversion score')

colormap gray
sgtitle([strrep(m.mouse,'_',' ') ' ' m.run])

% Image of mouse 
% webcampath = fullfile(strrep(m.CCDdir,'CCD','webcam'),m.run);
% frames = round([0:.25:.99]*m.movielength*60+1);
% mouse_frame = [];
% for i = 1:numel(frames)
%     try im1 = imresize(LoadFLIR(webcampath,frames(i),0),[1080 1440]);
%     catch; im1 = uint8(zeros(1080,1440)); end
%     try im2 = imresize(LoadFLIR(webcampath,frames(i),1),[1080 1440]);
%     catch; im2 = uint8(zeros(1080,1440)); end
%     mouse_frame(:,[1:1440]+1440*(i-1)) = [im1;im2];
% end
% subplot(sp_sz(1),sp_sz(2),sp_extent([5 5],[3 4],sp_sz))
% imagesc(mouse_frame); axis image; axis off
% title('mouse images at 0, 2.5, 5, 7.5 min')
