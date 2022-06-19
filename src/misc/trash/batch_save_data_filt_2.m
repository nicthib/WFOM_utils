%% preprocessing step
% this script loads, converts raw data and then saves jrgeco

%% Add these path
clear all; close all; clc;
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
addpath(genpath('/local_mount/space/dingus/1/Weihao/Resting_State/Utils'));
addpath /local_mount/space/dingus/1/Weihao/Resting_State/Code
save_root = '/local_mount/space/dingus/1/RS_analysis/jrgeco';
data_root = '/local_mount/space/dingus/1/cmdata/';

%% Pool
mouse = 'cm128';
days = ['3'];
%% Return path and stimlist [cell arrays organized by day]
for d = 1:length(days)   
    animaldayfolder = ['/local_mount/space/dingus/1/cmdata/' mouse '_' days(d)];
    [CCDpath{d},stimlist{d}]=getstimlist_RSanalysis(animaldayfolder);
end
%% one mask for all sessions
load(fullfile(save_root,mouse,'TF_w_wx.mat')); % registration file
load(fullfile(save_root,mouse,'1','mask_rs_v2_cm128_1.mat')) % bilateral
load(fullfile(save_root,mouse,'1','mask_rs_v1_cm128_1.mat')) % unilateral
%% load raw data 
bug_run = [];
for d = 1:length(days)
    
    day = days(d);
    mouseday = [mouse '_' day];    

    % go to each run in each day...
    for i_run = 1:numel(stimlist{d}) 
        run = stimlist{d}{i_run}(1:4);
        if ~strcmp(run, 'runA') % do every run except for runA  
            gopath = fullfile(data_root, mouseday,'CCD',run,stimlist{d}{i_run});
%             try
                %% Load data
                m = makem;
                m.dsf = 2;
                m.outputs = 'rgl';
                m.nrot = 2;
                m.corr_flicker = 1:3;
                m.loadpct = [0 1];
                m.BW = ones(256,256);
                % choose baseline:
                m.baseline = 60:1000;
                % load data in
                [m, data] = LoadData_v2(gopath,m);
                
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
                else
                    for i = 1:m.nLEDs
                        data.(m.LEDs{i}) = data.(m.LEDs{i})(:,:,1:end-5);% GET RID OF THE LAST 5 FRAMES (last spool)
                    end
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
                
                %%
                %bp = [0.1, 0.2, 0.3, 0.5, 0.8, 1, 3, 5]; % low pass portion of bp filter
                %lp = [0.02, 0.04, 0.06, 0.08] %0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28,0.3];
                %lp = [0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28,0.3, 1.1, 1.5, 2.1, 2.5, 3.1]
                lp = [0.1, 0.2, 0.35, 0.40, 0.45, 0.5, 0.75]; %this one next
                %bp = [1.1, 1.5, 2.1, 2.5, 3.1];
                %lp = 'unfiltered'
                for sf = 1: length(lp)
                    save_root = ['/local_mount/space/dingus/1/RS_analysis/jrgeco_filtered/',num2str(lp(sf))];
                    %save_root = ['/local_mount/space/dingus/1/RS_analysis/jrgeco_filtered/',lp];
                    disp(['save root =' num2str(lp(sf))])
                    
                    day = days(d);
                    mouseday = [mouse '_' day];
                    save_path = fullfile(save_root,mouse,day);
                    
                    % save to...
                    if ~exist(save_path)
                        mkdir(save_path);
                    end
                    
                    %% bp filter
%                     fs = m.framerate/3;
%                     f_n = [0, 0.02, bp(sf),bp(sf)+0.02]*2/fs; % stopband in normalized frame rate
%                     a = [1 0];
%                     b_hp = firpm(300,[0 f_n 1], [0 0 1 1 0 0]); % order is 300
%                     ss = size(data.jrgeco);
%                     non_nan_ind = reshape(~isnan(data.jrgeco(:,:,1)),[ss(1)*ss(2),1]);
%                     datahere = {'jrgeco'};
%                     for j = 1:numel(datahere)
%                         datareshape = (reshape(data.(datahere{j}),[ss(1)*ss(2),ss(3)]))';
%                         datareshape(:,non_nan_ind) = filtfilt(b_hp,a,datareshape(:,non_nan_ind));
%                         data.([datahere{j} '_filt']) = reshape(datareshape',[ss(1),ss(2),ss(3)]);
%                         clear datareshape
%                     end
%                     disp('done w hp')
                    %% LP filter
                    ss = size(data.jrgeco(:,:,1:end-15));
                    Rp=1; As=50;

                    fs = m.framerate/3;
                    f = [lp(sf),lp(sf)+0.02];

                    a = [1,0];
                    dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1),10^(-As/20)];
                    [M, f0, a0, weights] = firpmord(f, a, dev, fs);
                    b_lp = firpm(M+1,f0,a0, weights);
                    jrgeco_reshape = reshape(data.jrgeco(:,:,1:end-15),[ss(1)*ss(2),ss(3)]);
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
                    %data.jrgeco_filt = data.jrgeco;
    
                    
                    
                    %% now decompose jrgeco and save
                    ss= size(data.jrgeco);
                    datahere = {'chbt','jrgeco_filt'};
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
                    save(fullfile(save_path,stimlist{d}{i_run}),'C','S','expl','m','mask');
                    disp(['done saving ' stimlist{d}{i_run}])
                    %clear run gopath m C S expl
                    
                end
                
                

%             catch ME
%                 bug_run = [ bug_run; stimlist{d}{i_run}];
%             end
        end
    end
    
    
end
