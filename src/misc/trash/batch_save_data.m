%% preprocessing step
% this script loads, converts raw data and then saves jrgeco
%% 05/31/2020
% remove high-pass filtering part
% add new mask (Nic's most up-to-date)
% save root changed to behavior 2
%% 06/16/2020
% add cm128 support
% subplot(132);imagesc(data.green(:,:,100));
% subplot(131);imagesc(IDX.cm128.b1{1, 9});
% subplot(133);imagesc(mask)
% delete registration code
%% 06/19/2020
% create cm128 mask and IDX with midline sorted [ generated from
% make_IDX_cm125.m script 
% implement low-pass filter 0.02Hz order = 800
% register to cm128_1 frame by frame
%% 07/13/2020
% add cm126 registration
% dont want day 4 of cm126

%% 07/13/2020
% add cm126 registration
% dont want day 4 of cm126

%% 09/04/2020
% stop applying mask when loading data
% save output to chronos

dbstop if error
%% Add these path
clear all; close all; clc;
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
addpath(genpath('/local_mount/space/dingus/1/Weihao/Resting_State/Utils'));
addpath /local_mount/space/dingus/1/Weihao/Resting_State/Code
save_root = '/local_mount/space/chronos/1/WFOM_DATA_ANALYSIS/Weihao/analysis';
data_root = '/local_mount/space/dingus/1/cmdata/';

%% Pool
mouse = 'cm125';
days = ['1','2','3','5','6'];
%days = '7';
%% Return path and stimlist [cell arrays organized by day]
for d = 1:length(days)   
    animaldayfolder = ['/local_mount/space/dingus/1/cmdata/' mouse '_' days(d)];
    [CCDpath{d},stimlist{d}]=getstimlist_RSanalysis(animaldayfolder);
end


%% one mask for all sessions
% load(fullfile('/local_mount/space/dingus/1/RS_analysis/jrgeco','TF_w_wx.mat')); % registration file
% load(fullfile(save_root,mouse,['mask_' mouse '.mat'])) % unilateral
% load(fullfile(save_root,mouse,['IDX_' mouse '.mat'])) % midline sorted left and right
%load('/local_mount/space/dingus/1/RS_analysis/registration/TF_new.mat'); % load IDX struct

% 07/13/2020
% Now, all masks are the same - refers to register_cm126.m in jrgeco folder

load(fullfile('/local_mount/space/dingus/1/RS_analysis/behavior2','cm125',['mask_cm125.mat'])) % universal mask

registration file for cm126 
load(fullfile('/local_mount/space/dingus/1/RS_analysis/jrgeco/cm126','TF_w_wx_cm126.mat')); % within
load(fullfile('/local_mount/space/dingus/1/RS_analysis/jrgeco','TF_a_wx.mat')); % across

% registration file for cm125
% load(fullfile('/local_mount/space/dingus/1/RS_analysis/jrgeco/cm126','TF_w_wx_cm126.mat')); % within

%% load raw data 
bug_run = [];
for d = 1:length(days)
    
    day = days(d);
    mouseday = [mouse '_' day];    
    save_path = fullfile(save_root,mouse,day);
    
    % save to...
    if ~exist(save_path)
        mkdir(save_path);
    end
       
    % go to each run in each day...
    for i_run = 1:numel(stimlist{d}) 
        run = stimlist{d}{i_run}(1:4);
        if ~strcmp(run, 'runA') % do every run except for runA  
            gopath = fullfile(data_root, mouseday,'CCD',run,stimlist{d}{i_run});
            try
                %% Load data
                m = makem;
                m.dsf = 2;
                m.outputs = 'rgl';
                m.nrot = 2;
                %m.corr_flicker = 1:3;
                m.loadpct = [0 1];
                m.BW = ones(256,256);
                % choose baseline:
                m.baseline = 60:900;
                % load data in
                [m, data] = LoadData_v2(gopath,m);
                
                for i = 1:m.nLEDs
                    data.(m.LEDs{i}) = data.(m.LEDs{i})(:,:,1:end-5);% GET RID OF THE LAST 5 FRAMES (last spool)
                end
                %% registration
                if ~strcmp(day, '1') % register every day to day 1
                    ss= size(data.green);
                    m.BW = mask;
                    tic
                    for i = 1:m.nLEDs
                        disp(['registering ' m.LEDs{i}])
                        for iT = 1:ss(3)
                            temp(:,:,iT) = imwarp_fr(data.(m.LEDs{i})(:,:,iT),TF_w_wx.([m.mouse '_1']));
                        end
                        data.(m.LEDs{i}) = temp(:,:,1:end-5);% GET RID OF THE LAST 5 FRAMES (last spool)
                        data.(m.LEDs{i}) = data.(m.LEDs{i}).*repmat(imresize(m.BW,ss(1:2)),[1 1 size(data.(m.LEDs{i}),3)]);% apply mask
                    end
                    toc
                    clear temp
                    disp(['Done Registration ']);
                else % this is day 1
                    ss= size(data.green);
                    m.BW = mask;
                    tic
                    for i = 1:m.nLEDs
                        disp(['day 1 registering ' m.LEDs{i}])
                        for iT = 1:ss(3)
                            temp(:,:,iT) = imwarp_fr(data.(m.LEDs{i})(:,:,iT),TF_a_wx.([mouse 'tocm125']));
                        end
                        data.(m.LEDs{i}) = temp(:,:,1:end-5);% GET RID OF THE LAST 5 FRAMES (last spool)
                        data.(m.LEDs{i}) = data.(m.LEDs{i}).*repmat(imresize(m.BW,ss(1:2)),[1 1 size(data.(m.LEDs{i}),3)]);% apply mask
                    end
                    toc
                    clear temp
                    disp(['Done Registration ']);
                end
                %% PCA denoising
                m.PCAcomps = 500;
                ss= size(data.green);
                for i = 1:m.nLEDs
                    disp(['PCA ' m.LEDs{i}])
                    tmp = reshape(data.(m.LEDs{i}),[ss(1)*ss(2),ss(3)]);
                    m.nanidx = isfinite(tmp(:,1));% non nan indice
                    [C,S] = fsvd(tmp(m.nanidx,:),m.PCAcomps,2,1);
                    tmp(m.nanidx,:) = C*S;
                    data.(m.LEDs{i}) = reshape(tmp,ss);
                end
                clear tmp C S
                
                %% Convert Data
                disp('Converting Hemodynamics...')
                [data.chbo,data.chbr,data.chbt] = ...
                    convert_mariel_MIAO(mask.*data.green,mask.*data.red,'g','r',m.baseline,534);
                
                disp('Converting JRGECO...')
                Dg = 0.4;
                Dr = 1.115;
                ss = size(data.lime);
                data.jrgeco= repmat(mask,[1 1 ss(3)]).*data.lime./((abs(data.red-90).^Dr.*(abs(data.green.^Dg))));
                bgG=mean(data.jrgeco(:,:,m.baseline),3);
                data.jrgeco=data.jrgeco./repmat(bgG,[1 1 ss(3)])-1;
                %% hp filter
                fs = m.framerate/3;
                f_n = [0, 0.02]*2/fs; % stopband in normalized frame rate
                a = [1 0];
                b_hp = firpm(800,[0 f_n 1], [0 0 1 1]); % order is 800
                ss = size(data.jrgeco);
                non_nan_ind = reshape(~isnan(data.jrgeco(:,:,1)),[ss(1)*ss(2),1]);
                datahere = {'chbt','jrgeco'};
                for j = 1:numel(datahere)
                    datareshape = (reshape(data.(datahere{j}),[ss(1)*ss(2),ss(3)]))';
                    datareshape(:,non_nan_ind) = filtfilt(b_hp,a,datareshape(:,non_nan_ind));
                    data.([datahere{j} '_hp']) = reshape(datareshape',[ss(1),ss(2),ss(3)]);
                    clear datareshape
                end
                disp('done w hp')
% % %                 
                %% save sample movie
                vidfile = VideoWriter(fullfile(save_path,stimlist{d}{i_run}));
                vidfile.FrameRate = 20;
                open(vidfile)
                clear F
                close all;
                figure('Position',[0 0 620 300]);
                for i = 6000:6600
                    subplot(121)
                    imagesc(data.jrgeco_hp(:,:,i)); axis image off; colormap(gca,'gray');caxis([-1 1]*0.04);
                    title(num2str(i));
                    subplot(122)
                    imagesc(data.chbt_hp(:,:,i)); axis image off; title('cHbT');colormap(gca,'jet');caxis([-10 10]*1e-6);
                    F(i) = getframe(gcf);
                    pause(0.01);
                end
                writeVideo(vidfile,F(6002:6600))
                close(vidfile)
% %                 
                %% now decompose jrgeco and save
                ss= size(data.jrgeco);
                datahere = {'chbt','jrgeco','chbt_hp','jrgeco_hp'};
                clear C S expl
                for j = 1:numel(datahere)
                    disp(['PCA and saving ' datahere{j}])
                    datareshape = reshape(data.(datahere{j}),[ss(1)*ss(2),ss(3)]);
                    m.nanidx = isfinite(datareshape(:,1));% non nan indice
                    [C.(datahere{j}),S.(datahere{j}),expl.(datahere{j})] = fsvd(datareshape(m.nanidx,:),m.PCAcomps,2,1);
                end
                clear datareshape
                % save important things
                save(fullfile(save_path,stimlist{d}{i_run}),'C','S','expl','m','mask');

            clear run gopath m C S expl
            catch ME
                bug_run = [ bug_run; stimlist{d}{i_run}];
            end
        end
    end
    
    
end


