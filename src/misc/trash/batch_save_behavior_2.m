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
% getHfromkmeans first and then multi-bandpass data
% CROP OUT FIRST 5 AND LAST 15 FRAMES
%% update 06/16/2020
% load most up-to-date IDX map: add register code
% undo CROP OUT FIRST 5 AND LAST 15 FRAMES
%% update 6/20/2020
% now load standard cm128 mask
% add visualization of H_recon
%% update 6/28/2020
% load 42 comps map
% save to behavior 3 folder and H42_**
%% update 7/28/2020
% save chbt H too

%% Add these path
clear all; close all; clc;
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
addpath(genpath('/local_mount/space/dingus/1/Weihao/Resting_State/Utils'));
addpath /local_mount/space/dingus/1/Weihao/Resting_State/Code
save_root = '/local_mount/space/dingus/1/RS_analysis/behavior3';
data_root = '/local_mount/space/dingus/1/RS_analysis/behavior2';

%% Pool one mouse at a time
mouse = 'cm126';
days = ['1','2','3','5','6'];

%% Return path and stimlist [cell arrays organized by day]
clear animaldayfolder;
for d = 1:length(days)   
    day = days(d);
    animaldayfolder{d} = dir(fullfile(data_root,mouse,day,'*_stim_1.mat'));
end

%% get that one mask of this mouse and do kmeans - ONCE FOR ALL
 % day 1 run B

%load(fullfile(save_root,mouse,'sorted_idx.mat'),'IDX','IDX_sorted');
% load('/local_mount/space/dingus/1/RS_analysis/H_b_new/TF.mat','cm125');
% load('/local_mount/space/dingus/1/RS_analysis/registration/TF_new.mat')
% IDX = cm125.IDXs_b{1, 3};  %idx map w most comps
load(fullfile(data_root,mouse,['mask_' mouse '.mat'])) % unilateral
load(fullfile(data_root,mouse,['IDX_36comps_' mouse '.mat']),'IDX_sorted') % midline sorted left and right


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
        %         try
        filename = fullfile(mainfolder(i_run).folder,mainfolder(i_run).name);
        load(filename)
        % use 36 comps IDX
        load(fullfile(data_root,mouse,['IDX_36comps_' mouse '.mat']),'IDX_sorted') % midline sorted left and right

        % 256*256 -> 128*128
        T = size(S.jrgeco,2);
        jrgeco = nan(256*256,T);
        jrgeco(m.nanidx,:) = C.jrgeco*S.jrgeco;
        data.jrgeco = reshape(jrgeco,[256 256 T]);

        T = size(S.chbt,2);
        chbt = nan(256*256,T);
        chbt(m.nanidx,:) = C.chbt*S.chbt;
        data.chbt = reshape(chbt,[256 256 T]);
        %% CAUTION: CROP OUT FIRST 5 AND LAST 15 FRAMES HERE
        %         data.jrgeco(:,:,1:5)=[];
        %         data.jrgeco(:,:,end-14:end)=[];
        
        
        % % %         if mouse ~= 'cm125'
        % % %             IDX = imwarp_fr(IDX,TF_a.(['cm125to' mouse]));
        % % %         end
        % % %         if day ~= '1'
        % % %            IDX = imwarp_fr(IDX,TF_w.([mice{a} '_1_' days(b)]));
        % % %         end
        H = getHfromKmeans(data.jrgeco,IDX_sorted,0);
        H_chbt = getHfromKmeans(data.chbt,IDX_sorted,0);

        ss = size(data.jrgeco);
        ncomps = size(H,1);
        % reconstruct back to get average
        IDX_sorted_reshape = reshape(IDX_sorted,[ss(1)*ss(2),1]);
        H_recon = nan(ss(1)*ss(2),ss(3)); % pre-allocate
        H_chbt_recon = nan(ss(1)*ss(2),ss(3)); % pre-allocate
        for iN = 1:ncomps
            indice = IDX_sorted_reshape==iN;
            H_recon(indice,:) = repmat(H(iN,:),nnz(indice),1);
            H_chbt_recon(indice,:) = repmat(H_chbt(iN,:),nnz(indice),1);
        end
        H_recon = reshape(H_recon,[ss(1),ss(2),ss(3)]);
        H_chbt_recon = reshape(H_chbt_recon,[ss(1),ss(2),ss(3)]);
        try
            vidfile = VideoWriter(fullfile(save_path,['H_' mainfolder(i_run).name]));
            vidfile.FrameRate = 20;
            open(vidfile)
            clear F
            close all;
            % show comparison
            figure('Position',[0,0,1200,600]);
            for i =7000:7500
                subplot(131)
                imagesc(H_recon(:,:,i));axis image; colormap(gca,'gray'); caxis([-1 1]*0.05);
                title([num2str(i) ' H recon']);
                
                subplot(132)
                imagesc(data.jrgeco(:,:,i));axis image; colormap(gca,'gray'); caxis([-1 1]*0.05);
                title([num2str(i) ' jrgeco']);
                
                subplot(133)
                imagesc(H_chbt_recon(:,:,i));axis image; colormap(gca,'jet');caxis([-8 8]*1e-6);
                title([num2str(i) ' H chbt recon']);
                
                F(i) = getframe(gcf);
                pause(0.01);
            end
            writeVideo(vidfile,F(7002:7500))
            close(vidfile)
            
        catch
            disp(['BUGGGGGGGGG ',day ,' -',mainfolder(i_run).name])
            continue;
        end

        %%
        %         for i =10000:10980
        %         imagesc(data.jrgeco(:,:,i));colormap gray;caxis([-1 1]*0.05);title(num2str(i));pause(0.01);
        %         end
        
        %%
        
        % frequency parameters
        min_freq =  0.05;
        max_freq = 6;
        num_frex = 5;
        
        
        % frequencies vector
        frex = linspace(min_freq,max_freq,num_frex);
        frex = logspace(log10(min_freq),log10(max_freq),num_frex);
        
        clear spectrum
        % spectrum
        for i_f = 1:numel(frex)-1
            spectrum(i_f,:) = [frex(i_f) frex(i_f+1)];
        end
        
        % multi-bandpass filter
        pass_bw = 0.03;
        fs = m.framerate/3;
        
        % filtering
        clear b_bp H_bp
        for i = 1:size(spectrum,1)
            
            low = spectrum(i,1);
            high = spectrum(i,2);
            f_n = [low-pass_bw low high high+pass_bw]*2/fs;
            a = [1 0];
            b_bp(i,:) = firpm(500,[0 f_n 1],[0 0 1 1 0 0]); % order is 500
            
            temp = filtfilt(b_bp(i,:),a,H');
            H_bp(i).jrgeco = temp';
            
            
            clear low high f_n
        end
        clear temp
        clear C S expl jrgeco
        
        %%
        
        save(fullfile(save_path,['H_' mainfolder(i_run).name]),'H','H_chbt','H_bp','m','spectrum','b_bp','IDX_sorted')
        
        disp(['done saving day ',day ,' -',mainfolder(i_run).name])
        %         catch ME
        %             bug_run = [bug_run; ['day_' day '_' mainfolder(i_run).name]];
        %         end
    end
    
    
end