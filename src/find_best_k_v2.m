%% v1 draft - 04/06/2020 wx
%% v2 draft - 07/16/2021 nic
clear
load(fullfile('/local_mount/space/dingus/1/RS_analysis/preprocessing','params.mat'))
addpath(genpath('/local_mount/space/dingus/1/RS_analysis/code'))
idx = 1; data_all = [];
for j = 1%:numel(runs_for_analysis)
    for i = 1:5:numel(sessions_for_analysis)
        file_id = [sessions_for_analysis{i} '_' runs_for_analysis{j}];
        load_dir = findmousefolder(sessions_for_analysis{i},runs_for_analysis{j},'1');
        mouse_id = sessions_for_analysis{i}(1:5);
        m = m_default;
        m.IDX = m_all.(mouse_id).IDX_final; 
        m.refim = m_all.(mouse_id).refim;
        m.baselinewidth = 1200; 
        m.loadpct = [0 1];
        m.dofilter = 0;
        [m,data,~,~] = LoadData_RS(load_dir,m);
        m.BW = double(m.IDX>0); m.BW(m.BW==0) = NaN;
        data_all{idx} = data.jrgeco(:,:,m.baseline);
        data_all{idx} = data_all{idx}.*repmat(m.BW,[1 1 size(data_all{idx},3)]);
        data_all{idx} = reshape(data_all{idx},[256^2 size(data_all{idx},3)]);
        idx = idx+1;
    end
end

%% JRGECO
% range of different k
Ks = 2:6:150;
for i_data = 1:numel(data_all)
    clear M;
    for i_k = 1:numel(Ks)
        ncomps = Ks(i_k);
        disp(['now k=' num2str(ncomps)]);
        % do kmeans, max iteration 200 for convergence; repeat 10 times (u can adjust)
        [M(i_k).IDX , M(i_k).Centroid] = kmeans(data_all{i_data},ncomps,'MaxIter',150,'Distance','correlation','Replicates',1);
    end
    
    % WICD/ACCD
    % within-cluster distance (WICD)
    clear WICD ACCD
    for i_k = 1:numel(Ks)
        centroid = M(i_k).Centroid;
        idx = M(i_k).IDX;
        for k = 1:size(centroid,1) % for each cluster's centroid
            % distance b/t this centroid and every within-cluster vector
            distances = pdist2(centroid(k,:) , data_all{i_data}(idx==k,:),'correlation');
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
    
    % ratio
    ratio_jrgeco(i_data,:) = WICD./ACCD;
end

Ks = 2:3:150;
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
    
    % ratio
    ratio_chbt(i_data,:) = WICD./ACCD;
end


%%
Ks_jrgeco = 2:6:150;
%Ks_chbt = 2:3:150;
plot(Ks_jrgeco,ratio_jrgeco,'LineWidth',2);hold on;
%plot(Ks_chbt,ratio_chbt+0,'LineWidth',2);hold on;
ylabel('WICD/ACCD','FontSize',14)
xlabel('K','FontSize',14)
legend('run1-jrgeco','run2-jrgeco','run3-jrgeco','run1-chbt','run2-chbt','run3-chbt','FontSize',14)
set(gcf,'color','w');
