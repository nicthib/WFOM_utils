%% ML WORKFLOW
% OVERALL GOAL: predict movement using H_n
clear; close all
addpath /local_mount/space/revault/revault2/Ewoud/SRGAP2C/Colormap/
addpath(genpath('/local_mount/space/dingus/1/RS_analysis'))
load ERES_2.mat
nets = {1,10,[10 10 10]};
%% STEP 1: make model with mouse 1
for nn = 0:2
    expname = 'cm128_6';
    runs = 'BCD';
    filenames = {};
    for j = 1:numel(runs)
        filenames{end+1} = [expname '_run' runs(j) '_H.mat'];
    end
    mapidx = 6;
    H_zf = [];
    rotf_all = [];
    for i = 1:numel(filenames)
        load(filenames{i})
        H_z = zscore(cleanupH(H_n{mapidx}),[],2);
        H_zf = cat(2,H_zf,H_z);
        rotf_all = [rotf_all rotf];
    end
    H_zf = cleanupH(H_zf);
    rotf_all = zscore(rotf_all);
    % STEP 2: Train model
    ns = 3; % Number of timepoints
    sh = 0; % backwards shift
    netsize = [nets{nn+1}];
    net = feedforwardnet(netsize);
    net.trainParam.epochs = 1000;
    net.trainParam.max_fail = 12;
    [net, tr] = train(net,H_zf,rotf_all,'UseParallel','yes');
    disp(['Training for net ' mat2str(j) ' complete. tr = ' mat2str(numel(tr.epoch))])
    
    % STEP 3: Plot results for subject sufficient model
    t = linspace(0,numel(rotf_all)/20,numel(rotf_all));
    subplot(3,3,1+nn*3)
    pred =  net(H_zf);
    plot(t,rotf_all);
    hold on
    plot(t,pred);
    cc = corrcoef(rotf_all,pred);
    title(['R = ' mat2str(round(cc(2,1)*100)/100)])
    xlim([0 300]); ylim([-2 5]); xlabel('time(sec)')
    
    % STEP 4: Evaluate deficient model for mouse 2 and 3
    expname = {'cm125_6'};
    for e = 1:numel(expname)
        filenames = {};
        for j = 1:numel(runs)
            filenames{end+1} = [expname{e} '_run' runs(j) '_H.mat'];
        end
        H_zf = [];
        rotf_all = [];
        for i = 1:numel(filenames)
            load(filenames{i})
            H_z = zscore(cleanupH(H_n{mapidx}),[],2);
            H_zf = cat(2,H_zf,H_z);
            rotf_all = [rotf_all rotf];
        end
        rotf_all = zscore(rotf_all);
        subplot(3,3,e+1+nn*3)
        pred =  net(H_zf);
        plot(t,rotf_all);
        hold on
        plot(t,pred);
        cc = corrcoef(rotf_all,pred);
        title(['R = ' mat2str(round(cc(2,1)*100)/100)])
        xlim([0 300]); ylim([-2 5]); xlabel('time(sec)')
    end
end