%% v1 draft - 04/06/2020 wx
clear all; close all;
%% Load Example Data
% load(fullfile('/local_mount/space/dingus/1/Weihao/Resting_State/cmdata_CCD_analysis','cm128','RSP_cm128_3_runD.mat'),'C','expl','m');
% ncomps_svd = size(C.jrgeco,2);
% jrgeco_reduced = nan(256*256,ncomps_svd);
% jrgeco_reduced(m.nanidx,:) = C.jrgeco*diag(expl.jrgeco); % pixels by time

%% Load data
clear; load('params_checkpoint.mat')

%% JRGECO
sessions = {'cm124_2_runD','cm125_2_runD','cm126_2_runD','cm127_2_runD','cm128_2_runD'};
Ks = 5:10:300;
for i_data = 1:numel(sessions)
    mouse = sessions{i_data}(1:7);
    mouserun = sessions{i_data}(9:end);
    loadpath = findmousefolder(mouse,mouserun,1);
    m = m_defaults.preprocessing; % Default options
    m.outputs = 'n';
    [m,datatmp] = LoadData_v2(loadpath,m);
    
    BW = round(imwarp(m_source.BW,m_all.(mouse).TF,'OutputView',imref2d([m_source.sz m_source.sz])));
    BW(BW==0) = NaN;
    
    data_full = datatmp.jrgeco.*repmat(BW,[1 1 size(datatmp.jrgeco,3)]);
    data_full = reshape(data_full,[256*256 size(datatmp.jrgeco,3)]);
    
    clear M;
    [C,S,~,~,E] = pca(data_full);
    data_recon = S(:,1:1000)*diag(E(1:1000));
    for i_k = 1:numel(Ks)
        ncomps = Ks(i_k);
        disp(['now k=' num2str(ncomps)]);
        % do kmeans, max iteration 200 for convergence; repeat 10 times (u can adjust)
        [M(i_k).IDX , M(i_k).Centroid] = kmeans(data_recon,ncomps,'MaxIter',150,'Distance','correlation','Replicates',1);
    end
    
    % WICD/ACCD within-cluster distance (WICD)
    clear WICD ACCD
    for i_k = 1:numel(M)
        centroid = M(i_k).Centroid;
        idx = M(i_k).IDX;
        for k = 1:size(centroid,1) % for each cluster's centroid
            % distance b/t this centroid and every within-cluster vector
            distances = pdist2(centroid(k,:) , data_recon(idx==k,:),'correlation');
            % average them
            distancee(k) = mean(distances);
        end
        % average all clusters
        WICD(i_k) = mean(distancee);
    end
    
    % across-cluster distance
    clear distancee
    for i_k = 1:numel(M)
        centroid = M(i_k).Centroid;
        % mean paired distance across clusters
        ACCD(i_k) = mean(pdist(centroid,'correlation'));
    end
    
    % ratio
    ratio_jrgeco(i_data,:) = WICD./ACCD;
end

%% CHBT
for i_data = 1:numel(data)
    clear M;
    for i_k = 1:numel(Ks)
        ncomps = Ks(i_k);
        disp(['now k=' num2str(ncomps)]);
        % do kmeans, max iteration 200 for convergence; repeat 10 times (u can adjust)
        [M(i_k).IDX , M(i_k).Centroid] = kmeans(data(i_data).chbt_reduced,ncomps,'MaxIter',150,'Distance','correlation','Replicates',10);
    end
    
    %% WICD/ACCD
    % within-cluster distance (WICD)
    clear WICD ACCD
    for i_k = 1:numel(Ks)
        
        centroid = M(i_k).Centroid;
        idx = M(i_k).IDX;
        
        for k = 1:size(centroid,1) % for each cluster's centroid
            % distance b/t this centroid and every within-cluster vector
            distances = pdist2(centroid(k,:) , data(i_data).chbt_reduced(idx==k,:),'correlation');
            % average them
            distancee(k) = mean(distances);
        end
        % average all clusters
        WICD(i_k) = mean(distancee);
    end
    
    % across-cluster distance
    clear distancee
    for i_k = 1:numel(Ks)
        centroid = M(i_k).Centroid;
        % mean paired distance across clusters
        ACCD(i_k) = mean(pdist(centroid,'correlation'));
    end
    ratio_chbt(i_data,:) = WICD./ACCD;
end

%% Calculate variance explained, look for # > 95%
Cutoff = .95;
cla
PC = [];
for i = 1:size(ratio_jrgeco,1)
    Var = ratio_jrgeco(i,1:end-1) - ratio_jrgeco(i,2:end);% % calculate %variance explained
    PC(i,:) = cumsum(Var)/(ratio_jrgeco(1)-ratio_jrgeco(end));
    plot(Ks(2:end),PC(i,:),'LineWidth',1,'Color',[0 0 0 .3]); hold on;
end
plot([mean(bestK) mean(bestK)],[0 1],'k--')
plot(Ks,[mean(PC) 1],'LineWidth',3)
ylabel('% Variance Explained','FontSize',14)
xlabel('K','FontSize',14)
title(['Best K = ' mat2str(round(mean(bestK)))])
ylim([.5 1])