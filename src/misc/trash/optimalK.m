%% OPTIMAL K FOR RAW DATA
clear
root = '/local_mount/space/dingus/1/RS_analysis/cm_analysis';
allrunnames = getallexps(root,{'C','rotf'});
filename = [];
nruns = 20;

%%
[IDX,data.jrgeco] = sample_kmeans(allrunnames(cellfun(@(s)~isempty(regexp(s,'run[BCD]')),allrunnames)),nruns,11);
data.jrgeco = reshape(data.jrgeco,[256 256 size(data.jrgeco,ndims(data.jrgeco))]);

%% Perform kmeans
IDX_clk = kmeans(reshape(data.jrgeco(:,:,kmeans_idx),[sz^2,numel(kmeans_idx)]),size(seed_pts,1),'Distance','Correlation','Start',C);
IDX_clk = reshape(IDX_clk,[256 256]);
IDXs.source = IDX_clk;

%% Test against new data
% Load jrgeco
runname = allrunnames{randi(numel(allrunnames),1)}; runname = strrep(runname,'.mat','');
load([runname '.mat'],'m')
rotf = getrotf(m);
[data,H,~] = recon_SVD(runname,IDXs.source);
kmeans_idx = baselinefromrot(rotf,500,2); kmeans_idx = kmeans_idx(50:449);
IDXs.(runname) = kmeans(data.jrgeco(:,kmeans_idx),size(seed_pts,1),'Distance','Correlation','Start',H.jrgeco{1}(:,kmeans_idx));
IDXs.(runname) = reshape(IDXs.(runname),[256 256]);

%% Display
close all
titles = fieldnames(IDXs);
for i = 1:max(IDX2(:))
    for j = 1:numel(titles)
        subplot(2,3,j)
        IDXtmp = imgradient(IDXs.(titles{j})); IDXtmp(isnan(IDXtmp)) = 0; IDXtmp(IDXtmp~=0) = 1;
        comptmp = IDXs.(titles{j})==i;
        imagesc(comptmp + IDXtmp); hold on
        %rectangle('Position',[seed_pts(i,1)-bx seed_pts(i,2)-bx bx*2 bx*2],'Facecolor','w','EdgeColor','None')
        axis image; axis off; title(strrep(titles{j},'_',' '))
    end
    pause
end





