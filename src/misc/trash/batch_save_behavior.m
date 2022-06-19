%% v1  - 04/30/2020 wx
% batch convert jrgeco to
% (1) bandpass H
% % % % % % % % % % % % (2) bandpass raw data [too large to save]
% (3) motion energy PC in differnt parts of body
% (4) m data
% (5) camera data (in different file)

%% CAUTION
%% CAUTION
%% CAUTION: CROP OUT FIRST AND LAST 500 FRAMES!!!
%% SO PLEASE BE CAREFUL IN THE BEHAVIOR DATA [1, 500] [end-500,end]

%% update 05/31/2020
% load most up-to-date IDX map

%% Add these path
clear all; close all; clc;
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
addpath(genpath('/local_mount/space/dingus/1/Weihao/Resting_State/Utils'));
addpath /local_mount/space/dingus/1/Weihao/Resting_State/Code
save_root = '/local_mount/space/dingus/1/RS_analysis/behavior2';
data_root = '/local_mount/space/dingus/1/RS_analysis/jrgeco';

%% Pool one mouse at a time
mouse = 'cm125';
days = ['1','2','3','4','5','6'];
%% Return path and stimlist [cell arrays organized by day]
clear animaldayfolder;
for d = 1:length(days)   
    day = days(d);
    animaldayfolder{d} = dir(fullfile(data_root,mouse,day,'*_stim_1.mat'));
end

%% get that one mask of this mouse and do kmeans - ONCE FOR ALL
 % day 1 run B
% % % % % main_folder = animaldayfolder{1};
% % % % % filename = fullfile(mainfolder(1).folder,mainfolder(1).name);
% % % % % load(filename);
% % % % % % do kmeans 
% % % % % T = size(S.jrgeco,2);
% % % % % jrgeco = nan(256*256,T);
% % % % % jrgeco(m.nanidx,:) = C.jrgeco*S.jrgeco;
% % % % % jrgeco = reshape(jrgeco,[256 256 T]);
% % % % % jrgeco = squeeze(nanmean(nanmean(reshape(jrgeco,[2,128,2,128,T]),3),1)); 
% % % % % jrgeco_reshape = reshape(jrgeco,[128*128, T]);
% % % % % m.nanidx_ds = ~isnan(jrgeco_reshape(:,1));
% % % % % jrgeco_reshape = jrgeco_reshape(m.nanidx_ds,:);
% % % % % jrgeco_reduced = fsvd(jrgeco_reshape,500,2,true);
% % % % % 
% % % % % ncomps = 100;
% % % % % IDX  = kmeans(jrgeco_reduced,ncomps,'MaxIter',200,'Distance','correlation','Replicates',200);
% % % % % real_idx = zeros(128*128,1);
% % % % % real_idx(m.nanidx_ds)=IDX;
% % % % % real_idx = reshape(real_idx,[128 128]);
% % % % % figure()
% % % % % imagesc(real_idx);
% % % % % 
% % % % % IDX=real_idx;
% % % % % [IDX_sorted,I0] = kmeans_sorted_shk_wx(ncomps,real_idx);
% % % % % imagesc(IDX_sorted);
% % % % % showIDX(IDX_sorted)
% % % % % colormap([1,1,1;'jet']);
% % % % % if ~exist(fullfile(save_root,mouse,'1'))
% % % % %     mkdir(fullfile(save_root,mouse,'1'));
% % % % % end
% % % % % save(fullfile(save_root,mouse,'sorted_idx'),'IDX','IDX_sorted');
% % % % % saveas(gcf,fullfile(save_root,mouse,'kmeans_idx.jpg'));

%load(fullfile(save_root,mouse,'sorted_idx.mat'),'IDX','IDX_sorted');
load('/local_mount/space/dingus/1/RS_analysis/H_b_new/TF.mat','cm125');
IDX_sorted = cm125.IDXs_b{1, 3};
%% loading data in
bug_run = [];
for d = 1:length(days)
    day = days(d);
    mainfolder = animaldayfolder{d};
    
    save_path = fullfile(save_root,mouse,day);
    
    % save to...
    if ~exist(save_path)
        mkdir(save_path);
    end
    
    %% get runs
    for i_run = 1:numel(mainfolder)
        try
        filename = fullfile(mainfolder(i_run).folder,mainfolder(i_run).name);
        load(filename)

        % 256*256 -> 128*128
        T = size(S.jrgeco,2);
        jrgeco = nan(256*256,T);
        jrgeco(m.nanidx,:) = C.jrgeco*S.jrgeco;
        data.jrgeco = reshape(jrgeco,[256 256 T]);

%         % downsample to 128*128
%         T = size(S.jrgeco,2);
%         data.jrgeco = squeeze(nanmean(nanmean(reshape(jrgeco,[2,128,2,128,T]),3),1));
        
        %%
        for i =10000:10980
        imagesc(data.jrgeco(:,:,i));colormap gray;caxis([-1 1]*0.05);title(num2str(i));pause(0.01);
        end
        %% CAUTION: CROP OUT FIRST AND LAST 500 FRAMES HERE
        %data.jrgeco(:,:,1:500)=[];
        data.jrgeco(:,:,end-499:end)=[];
        
        %%
        
        % frequency parameters
        min_freq =  0.05;
        max_freq = 9;
        num_frex = 8;
        
        
        % frequencies vector
        frex = linspace(min_freq,max_freq,num_frex);
        frex = logspace(log10(min_freq),log10(max_freq),num_frex);
        
        clear spectrum
        % spectrum
        for i_f = 1:numel(frex)-1
            spectrum(i_f,:) = [frex(i_f) frex(i_f+1)];
        end
        
        % multi-bandpass filter
        pass_bw = 0.02;
        fs = m.framerate/3;
        clear b_bp Mdata
        % pre-filtering
        ss = size(data.jrgeco);
        non_nan_ind = reshape(~isnan(data.jrgeco(:,:,1)),[ss(1)*ss(2),1]);
        datareshape = (reshape(data.jrgeco,[ss(1)*ss(2),ss(3)]))'; % shape ready for filtfilt
        temp = nan(size(datareshape));
        % filtering
        clear b_bp Mdata
        for i = 1:size(spectrum,1)
            
            low = spectrum(i,1);
            high = spectrum(i,2);
            f_n = [low-pass_bw low high high+pass_bw]*2/fs;
            a = [1 0];
            b_bp(i,:) = firpm(350,[0 f_n 1],[0 0 1 1 0 0]); % order is 350
            
            temp(:,non_nan_ind) = filtfilt(b_bp(i,:),a,datareshape(:,non_nan_ind));
            Mdata(i).jrgeco = reshape(temp',[ss(1),ss(2),ss(3)]);
            
            
            clear low high f_n
        end
        clear temp datareshape
        clear C S expl jrgeco
        
        %% dim reduction to extract kmeans regional wavelet signal and phase
        
        data_pool = {'jrgeco'};
        clear H
        for i =1:numel(data_pool)
            for i_freq = 1:numel(Mdata)
                H(i_freq).(data_pool{i}) = getHfromKmeans(Mdata(i_freq).(data_pool{i}),IDX_sorted,0);
            end
        end
        clear datahere
        %%

        save(fullfile(save_path,mainfolder(i_run).name),'H','m','spectrum','b_bp','IDX_sorted')

        disp(['done saving day ',day ,' -',mainfolder(i_run).name])
        catch ME
            bug_run = [bug_run; ['day_' day '_' mainfolder(i_run).name]];
        end
    end
    
    
end