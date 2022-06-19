%% 07/03/2020
% load jrgeco data in behavior2 and lp filter them in 15 different
% frequency
% saving to jrgeco_filtered folder - in different spatial comps
%% Add these path
clear all; close all; clc;
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
addpath(genpath('/local_mount/space/dingus/1/Weihao/Resting_State/Utils'));
addpath /local_mount/space/dingus/1/Weihao/Resting_State/Code
save_root = '/local_mount/space/dingus/1/RS_analysis/jrgeco_filtered/36comps';%%%%%%%%%%%%%%%%%%%%%%% REMEMBER TO CHANGE 
data_root = '/local_mount/space/dingus/1/RS_analysis/behavior2';

%% Pool one mouse at a time
mouse = 'cm125';
days = ['1','2','3','4','5','6'];
%days = ['1'];
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
load(fullfile(data_root,mouse,['IDX_36comps_' mouse '.mat'])) % midline sorted left and right

%% spectrum config
min_freq =  0.1;
max_freq = 9;
num_frex = 15;
passband_width = 0.03;

% frequencies vector
frex = linspace(min_freq,max_freq,num_frex);
frex = logspace(log10(min_freq),log10(max_freq),num_frex);

clear spectrum spectrum_string
% spectrum
for i_f = 1:numel(frex)
    spectrum(i_f,:) = [frex(i_f) frex(i_f)+passband_width];
    spectrum_string{i_f} = [num2str(frex(i_f)) '_' num2str(frex(i_f)+passband_width)];
end

%% loading data in
bug_run = [];
for d = 1:length(days)
    day = days(d);
    mainfolder = animaldayfolder{d};
    
    
    %% get runs
    for i_run = 1:numel(mainfolder)
%         try
        filename = fullfile(mainfolder(i_run).folder,mainfolder(i_run).name);
        load(filename)
        
        % 256*256 -> 128*128
        T = size(S.jrgeco,2);
        jrgeco = nan(256*256,T);
        jrgeco(m.nanidx,:) = C.jrgeco*S.jrgeco;
        data.jrgeco = reshape(jrgeco,[256 256 T]);
        %% CAUTION: CROP OUT FIRST 5 AND LAST 15 FRAMES HERE
        %         data.jrgeco(:,:,1:5)=[];
        %         data.jrgeco(:,:,end-14:end)=[];
        
        for i_f = 1: length(spectrum)
            savepath_root = ['/local_mount/space/dingus/1/RS_analysis/jrgeco_filtered/',spectrum_string{i_f}];
            %save_root = ['/local_mount/space/dingus/1/RS_analysis/jrgeco_filtered/',lp];
            disp(['save root =' spectrum_string{i_f}])
            
            day = days(d);
            mouseday = [mouse '_' day];
            save_path = fullfile(savepath_root,mouse,day);
            
            % save to...
            if ~exist(save_path)
                mkdir(save_path);
            end
            
            %% LP filter
            ss = size(data.jrgeco);
            Rp=1; As=50;
            
            fs = m.framerate/3;
            f = spectrum(i_f,:);
            
            a = [1,0];
            dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1),10^(-As/20)];
            [M, f0, a0, weights] = firpmord(f, a, dev, fs);
            b_lp = firpm(M+1,f0,a0, weights);
            jrgeco_reshape = reshape(data.jrgeco,[ss(1)*ss(2),ss(3)]);
            ss = size(data.jrgeco);
            non_nan_ind = reshape(~isnan(data.jrgeco(:,:,1)),[ss(1)*ss(2),1]);
            datahere = {'jrgeco'};
            for j = 1:numel(datahere)
                datareshape = (reshape(data.(datahere{j}),[ss(1)*ss(2),ss(3)]))';
                datareshape(:,non_nan_ind) = filtfilt(b_lp,a,datareshape(:,non_nan_ind));
                data.([datahere{j} '_filt'])  = reshape(datareshape',[ss(1),ss(2),ss(3)]);
                clear datareshape
            end
            disp('done w lp')
            
            %% CROP FIRST AND LAST 500 FRAMES
            data.jrgeco_filt(:,:,1:500)=[];
            data.jrgeco_filt(:,:,end-499:end)=[];
            
            vidfile = VideoWriter(fullfile(save_path,mainfolder(i_run).name));
            vidfile.FrameRate = 20;
            open(vidfile)
            clear F
            close all;
            figure('Position',[0 0 620 300]);
            for i = 6000:6600
                subplot(121)
                imagesc(data.jrgeco(:,:,i)); axis image off; colormap(gca,'gray');caxis([-1 1]*0.04);colorbar;
                title(num2str(i));
                subplot(122)
                imagesc(data.jrgeco_filt(:,:,i)); axis image off; colormap(gca,'gray');caxis([-1 1]*0.02);colorbar;
                F(i) = getframe(gcf);
                pause(0.01);
            end
            writeVideo(vidfile,F(6002:6600))
            close(vidfile)
            %% now decompose jrgeco and save
            ss= size(data.jrgeco_filt);
            datahere = {'jrgeco_filt'};
            clear C S expl
            m.PCAcomps = 500;
            for j = 1:numel(datahere)
                disp(['PCA and saving ' datahere{j}])
                datareshape = reshape(data.(datahere{j}),[ss(1)*ss(2),ss(3)]);
                m.nanidx = isfinite(datareshape(:,1));% non nan indice
                [C.(datahere{j}),S.(datahere{j}),expl.(datahere{j})] = fsvd(datareshape(m.nanidx,:),m.PCAcomps,2,1);
            end
            clear datareshape
            % save important things
            save(fullfile(save_path,mainfolder(i_run).name),'C','S','expl','m','mask');
            disp(['done saving ' mainfolder(i_run).name])
            %clear run gopath m C S expl
            
        end
        
                
%         catch ME
%             bug_run = [bug_run; ['day_' day '_' mainfolder(i_run).name]];
%             disp(['BUG ' bug_run]);
%         end
    end
    
    
end
