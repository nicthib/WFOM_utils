%% Loading files, functions and Data
clear
load('/local_mount/space/dingus/1/RS_analysis/preprocessing/params.mat')

%% Creating Global Functional Connectivity Map 
%close all
H_all = [];
for n = 1:numel(runnames_RS)
        fn = runnames_RS{n};
        load(fullfile(H_dir,fn)) % Loads the run
        rotf = smooth(m.rotf,100); % Smoothing the running data in the run for some reason
        H_all = [H_all H.jrgeco(:,rotf < .1)]; % per mouse, aggregates all the neural data in each ROI if the smoothed running is above 0.1
end
FC = corr(H_all'); 
figure
[ord,IDX_sub,FCIDX_sub] = FC_cluster(H_all,m.IDX,6); % Sorts the correlation matrix and puts the H data into the shuffled ROI map from earlier
% add networks variable
% Add save
% networks --> FC_groups

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
cmap = lines(7); cmap = cmap(2:7,:); % To change the colormap:
        % camp = colormap("pickColorMap"("# of different colors you want"))
subplot(223)
IDX_reord = zeros(size(IDXB));
for i = 1:max(IDXB(:))
   IDX_reord(IDXB == ord(i)) = i; 
end
subregions = histcounts(FCIDX_sub)/2; 
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
M2_regions = cumsum([0 histcounts(FCIDX_sub)/2]);% Not sure what histcounts is, but IDs largest ROI in each cluster after shuffling; why not handselect?
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