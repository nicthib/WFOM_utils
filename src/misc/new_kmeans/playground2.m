%% OPTIMAL K FOR RAW DATA
clear; root = '/local_mount/space/dingus/1/RS_analysis/cm_analysis';
allrunnames = getallexps(root,{'C'}); allrunnames = allrunnames';
runnames = allrunnames(cellfun(@(s)~isempty(regexp(s,'cm12[45678]_[345678]_run[BCDHIJ]')),allrunnames));
nruns = 20;
K = 80;
[IDX,jrgeco] = sample_kmeans(runnames,nruns,K);
IDX = reshape(IDX,[256 256]);

%%
for i = 1:numel(runnames)
    [jrgeco,m] = SVDtoGECO(runnames{i},1:500);
    if ischar(m.bl)
        IDX{i} = [];
    else
        kmeans_idx = baselinefromrot(m.rot,500,1); kmeans_idx = kmeans_idx(50:449);
        IDX{i} = reshape(kmeans(jrgeco(:,kmeans_idx),K,'Distance','Correlation'),[256 256]);
    end
    i
end

%% Creating decomposed brainmap with 16 components
IDX(isnan(IDX)) = 0;
[L,R,C,IDXout] = splitIDX(IDX);
for i = 1:size(IDXout,1)
    for j = 1:size(IDXout,2)
        if IDXout(i,j) == 0
            IDXout(i,j) = mode(mode(IDXout(i-3:i+3,j-3:j+3)));
        end
    end
end
[L,R,C,IDXout2] = splitIDX(IDXout);
compnames = {'L1','L2','L3','L4','L5','L6','L7',...
             'C1','C2',...
             'R1','R2','R3','R4','R5','R6','R7'};
for i = 1:max(IDXout2(:))
    IDXtmp = IDXout2==i;
    IDXtmp = double(IDXtmp);
    IDXtmp(IDXtmp==0) = NaN;
    IDXL1.(compnames{i}) = IDXtmp;
end

%% Show
spn = [13 14 9 10 5 6 1 2 3 16 15 12 11 8 7 4];
fn = fieldnames(IDXL1);
IDXg = imgradient(IDXout2); IDXg(isnan(IDXg)) = 0;
for i = 1:numel(fn)
    subplot(4,4,spn(i))
    imshowpair(IDXL1.(fn{i}),logical(IDXg)); axis image; axis off
    title(fn{i})
end

%% Level 2
K_L2 = 5;
for i = 1:numel(fn)
    IDX_test = double(IDXL1.(fn{i}));
    IDX_test(IDX_test==0) = NaN;
    clust_data = data.jrgeco.*repmat(IDX_test,[1 1 size(data.jrgeco,3)]);
    clust_data = reshape(clust_data,[256^2,size(data.jrgeco,3)]);
    IDXL2.(fn{i}) = kmeans(clust_data,K_L2,'Distance','Correlation');
    IDXL2.(fn{i}) = reshape(IDXL2.(fn{i}),[256 256]);
end

%% Show
close all
for i = 1:7
    c = num2str(i); 
    subplot(2,7,i)
    imagesc(IDXL2.(['L' c])); axis image; axis off; caxis([0 6])
    title(['L' c])
    subplot(2,7,i+7)
    imagesc(fliplr(IDXL2.(['R' c]))); axis image; axis off; caxis([0 6])
    title(['R' c'])
end


%% Test against new data
runname = allrunnames{randi(numel(allrunnames),1)}; runname = strrep(runname,'.mat','');
%runname='cm125_3_runD';
load([runname '.mat'],'m')
rotf = getrotf(m);
[data,H,~] = recon_SVD(runname,IDXL2);
data.jrgeco = reshape(data.jrgeco,[256 256 size(data.jrgeco,ndims(data.jrgeco))]);
data.jrgeco = data.jrgeco(:,:,1:end-10); rotf = rotf(1:end-10);
kmeans_idx = baselinefromrot(rotf,500,2); kmeans_idx = kmeans_idx(50:449);

for i = 1:numel(fn)
    clust_data = data.jrgeco(:,:,kmeans_idx).*repmat(IDXL1.(fn{i}),[1 1 numel(kmeans_idx)]);
    clust_data = reshape(clust_data,[256^2 numel(kmeans_idx)]);
    IDXL3.(fn{i}) = kmeans(clust_data,K_L2,'Distance','Correlation','Start',H.jrgeco.(fn{i})(:,kmeans_idx));
    IDXL3.(fn{i}) = reshape(IDXL3.(fn{i}),[256 256]);
end
IDXL4.(runname) = IDXL3;

%% Show
figure
%run1 = IDXL4.cm125_3_runD;
%run2 = IDXL4.cm127_2_runB;
IDXtest = IDXL4.cm128_2_runF;
for i = 1:7
    c = num2str(i); 
    subplot(2,7,i)
    imagesc(IDXtest.(['L' c])); axis image; axis off; caxis([0 K_L2])
    %tmp = run1.(['L' c])-run2.(['L' c]); tmp(isnan(tmp)) = 0;
    %imagesc(logical(tmp)); axis image; axis off; caxis([0 K_L2])
    title(['L' c])
    subplot(2,7,i+7)
    imagesc(fliplr(IDXtest.(['R' c]))); axis image; axis off; caxis([0 K_L2])
    %tmp = run1.(['R' c])-run2.(['R' c]); tmp(isnan(tmp)) = 0;
    %imagesc(logical(tmp)); axis image; axis off; caxis([0 K_L2])
    title(['R' c'])
end



