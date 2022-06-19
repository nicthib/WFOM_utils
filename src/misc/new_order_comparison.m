%% Loading files, functions and Data
clear
Hdir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; % The base folder we are pulling data from
addpath(genpath('/local_mount/space/dingus/1/RS_analysis/code'));
load('/local_mount/space/dingus/1/RS_analysis/Draft/params.mat')

% Creates cell with each element = list of mouse runs. (sep by mouse)
for i = 1:4
    runnames_singlemouse{i} = runnames_RS(cellfun(@(s) ~isempty(regexp(s,['cm12' mat2str(i+4)])),runnames_RS)); % Parse out runs that are 'mouse'
end

load(fullfile(Hdir,runnames_singlemouse{1}{1})) % load mouse data cm125_6_runJ

% Creating kmeans clustered ROI map
IDX_B = m.IDX; % Bilateral subregion map
IDX_B(m.IDX > 46) = IDX_B(m.IDX > 46)-46; % Combines pre-labelled ROI regions across central sulcus (92 --> 46 ROIs)
[IDX_B,I] = sortkmeans(IDX_B); % Creates an ROI based on a predetermined kmeans clustered correlation matrix
IDX_U = IDX_B; % Unilateral subregion map
IDX_B(m.IDX > 46) = IDX_B(m.IDX > 46)+46; % Splits the shuffled map into 92 ROIs
HI = [I;I+46]; % in order to find regions, both IDX and H must match. So we reordered IDX, we need to also reorder H later.

%% Creating Global Functional Connectivity Map 
close all
H_all = [];
FC_sr = []; FC_m = [];
for n = 1:4
    H_m = [];
    for r = 1:numel(runnames_singlemouse{n})
        fn = runnames_singlemouse{n}{r};
        load(fullfile(Hdir,fn)) % Loads the run
        rotf = smooth(m.rotf,100); % Smoothing the running data in the run for some reason
        H_m = [H_m H.jrgeco(:,rotf < .1)]; % per mouse, aggregates all the neural data in each ROI if the smoothed running is less than 0.1
        FC_sr{n}(:,r) = reshape(triu(corr(H.jrgeco(:,rotf < .1)')),[92^2 1]); % Mouse level FC matrix
    end
    FC_m(:,:,n) = corr(H_m'); % Mouse level FC matrix
    H_all = [H_all H_m]; % Combines all the neural data (why not done in inner loop?)
end

H_all_old_ord = H_all;
H_all_new_ord = H_all(HI,:);
FC_old_ord = corr(H_all_old_ord'); % Correlating all the neural data......make an average correlation map? I guess it's better than making one for each mouse? but why?
FC_new_ord = corr(H_all_new_ord'); %

[ord,IDX_sub,FC_IDX_sub] = sortH(H_all_new_ord,IDX_U,6); % Sorts the correlation matrix and puts the H data into the shuffled ROI map from earlier

%% Show trial specific FC corrs
FC_simscore = mean(corr(cell2mat(FC_sr)).^2);
score_mat = nan(4,47);
for n = 1:4
    subplot(1,4,n)
    idx = 1:numel(runnames_singlemouse{n});
    score_mat(n,idx) = FC_simscore(idx);
    FC_simscore(idx) = [];
end
subplot(211)
imagesc(score_mat'); axis image; axis off; caxis([.5 1]); colormap jet; colorbar
title('FC similarity (r squared) distribution')

subplot(212)
errorbar([1:4],nanmean(score_mat,2),nanstd(score_mat,[],2))

%% Reorder IDX: With Pre and Post Shuffling ROI maps;
figure
load(fullfile(Hdir,runnames_singlemouse{1}{end})) % Call cm125_6_runJ

% Shows old IDX map
cmap = lines(7); cmap = cmap(2:6,:);
close all
subplot(221)
subregions = [9 16 9 5 7];
showIDX_labels_subregions(m.IDX,subregions,cmap); % shows old clustering
ticklabels = {'M1/M2 L','SS L','S1/S2 L','RS L','Vis L',...
              'M1/M2 R','SS R','S1/S2 R','RS R','Vis R'};
title('Old order - region map')

% Showing FC (pre-shuffled)
subplot(222)
show_state_centroid(FC,ticklabels,subregions,[cmap;cmap])
colorbar
title('Old order - FC map')

% Reordered IDX map
cmap = lines(7); cmap = cmap(2:7,:);
subplot(223)
IDX_reord = zeros(size(IDX_B));
for i = 1:max(IDX_B(:))
   IDX_reord(IDX_B == ord(i)) = i; 
end
subregions = histcounts(FC_IDX_sub)/2; 
showIDX_labels_subregions(IDX_reord,subregions,cmap); %Show new clustering
ticklabels = {'Vis/RS L','M1/M2/SS L','SS/barrel L','SS/aud L','SS L','M2/M1 L',...
              'Vis/RS R','M1/M2/SS R','SS/barrel R','SS/aud R','SS R','M2/M1 R'};
ticklabels = ticklabels([[6 2 5 3 4 1], [6 2 5 3 4 1]+6]);
title('New order - region map')

% Showing FC after shuffling
subplot(224)
show_state_centroid(FC(HI(ord),HI(ord)),ticklabels,subregions,[cmap;cmap])
colorbar
title('New order - FC map')

% Calculating the average correlation
final_ord = HI(ord); % What is this? And why isn't it used?

M1 = FC;
M1_regions = cumsum([0 9 16 9 5 7]); % Largest ROI in Each Cluster (before shuffling); 
    % why not handselect the regions? Why use cumulative sum
M1_corr = [];

M2 = FC(HI(ord),HI(ord)); % A sorted correlation map with new order (but usng all the neural data)
M2_regions = cumsum([0 histcounts(FC_IDX_sub)/2]);% Not sure what histcounts is, but IDs largest ROI in each cluster after shuffling; why not handselect?
M2_corr = [];

%Isolating the Correlation in each cluster
for i = 1:numel(M1_regions)-1
    idx = M1_regions(i)+1:M1_regions(i+1); %Creates vector of ROIs in every cluster
    tmpcorr = M1(idx,idx); % isolates correlation matrix of the cluster
    M1_corr = [M1_corr; tmpcorr(:)]; % Vectorizes the lower triangle of isolated correlation matrix
end

for i = 1:numel(M1_regions)-1
    idx = M2_regions(i)+1:M2_regions(i+1);
    tmpcorr = M2(idx,idx);
    M2_corr = [M2_corr; tmpcorr(:)];
end

%Calculating mean
disp(['Within corr original = ' mat2str(round(mean(M1_corr),2))])
disp(['Within corr new = ' mat2str(round(mean(M2_corr),2))])
%save('new_order.mat','final_ord','subregions','ticklabels','cmap')