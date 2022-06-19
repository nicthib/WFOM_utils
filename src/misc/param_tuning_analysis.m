%% tune parameters -wx
% what params to tune in this script?
% -window size 
% -step size
% -window type

% p.s. tune number of clusters using gap statistics
% https://web.stanford.edu/~hastie/Papers/gap.pdf

%% Correlation analysis
% Setup, options
clear all; close all; clc; 
addpath(genpath('/local_mount/space/dingus/1/Weihao/Resting_State/Behavior/Code/Utils')) 
addpath(genpath('/local_mount/space/dingus/1/RS_analysis/code')) 


%%

opts.mapidx = 17;

exp_root = '/local_mount/space/dingus/1/RS_analysis/H_b';
expname = {'cm125_5','cm125_6'};%,'cm126_5','cm126_6','cm128_4','cm128_6'};
runs = 'BCD';
runnames = {};
for i = 1:numel(expname)
    for j = 1:numel(runs)
        runnames{end+1} = fullfile(exp_root,expname{i}(1:5),expname{i}(end),[expname{i} '_run' runs(j) '_H.mat']);
    end
end

% load H
H = []; delrows = []; % Initialize H and delrows
disp(sprintf(['Getting correlation states for \n' char(join(runnames,'\n'))]))
if isfield(opts,'H_ord')
    H_ord = opts.H_ord;
else
    load(runnames{1});
    H_ord = 1:size(H_n{opts.mapidx},1);
end
for i = 1:numel(runnames)
    load(runnames{i}) % Load data
    [~,delrows] = cleanupH(H_n{opts.mapidx}(H_ord,1:end-5),delrows); % Full H
end
delrows = unique(delrows);
for i = 1:numel(runnames)
    load(runnames{i}) % Load data
    [Htmp, ~] = cleanupH(H_n{opts.mapidx}(H_ord,1:end-5),delrows); % Full H
    H{end+1} = padarray(Htmp,[0 5],'pre');
end


%% GRID SEARCH - hyperparam tuning

% // 3D grid search so that i can visualize it using 3d surface\\
win_sizes = [35:20:435];
step_sizes = [5:20:105];
is_taper_window = [1 0];

opts.mapidx = 17; 
% mapidx 6: 26 comps; 
% mapidx 16: 76 comps; 
% mapidx 17: 81 comps;
% mapidx 20: 96 comps;
opts.nreps = 50; 
opts.sig = 0; % square if sig = 0, tapered otherwise 
opts.dokmeans = 1; 
tic
clear s c d1 dist_2_cent WICD
for i = 1:numel(win_sizes)
    for j = 1:numel(step_sizes)
        for k = 1:numel(is_taper_window)
            fprintf('now: win_sizes %d step_sizes %d is_taper_window %d \n',...
                win_sizes(i) , step_sizes(j) , is_taper_window(k)) 
            
            opts.ww = win_sizes(i);
            opts.ds = step_sizes(j);%floor(opts.ww/2); is this step size?
            opts.k = num_ks(k); 
            [~,c,~,dist_2_cent] = getcorrstates_wx(H,opts);
            %% compute WICD
            temp=cell2mat(dist_2_cent);
            WICD(i,j,k)=mean(temp(:));
            %% compute ACCD
            temp = pdist(c);
            ACCD(i,j,k)=mean(temp(:));
            clear temp
        end
    end   
end
toc
%% compute WICD


%% link orders (for long vs short)
[ordc, ~] = orderstates({c1{1},c2{1}},1:k);
for i = 1:numel(c2)
    c2{i} = c2{i}(ordc(2,:),:);
end
% combine
c = {c1{:},c2{:}};

%% Show corr maps 
[ord, distmin] = orderstates(c,1:k);
f1 = figure;% f2 = figure;
for i = 1:size(c,2)
    for j = 1:k
        figure(f1)
        subplot(size(c,2),k,(i-1)*k+ord(i,j))
        imagesc(reshape(c{i}(j,:),[r r]))
        caxis([0 1]); 
        axis image; axis off
        title(['S' num2str(j), ', dist = ' mat2str(round(distmin(i,j)/r,2))])
        colormap(jet)
        %colorbar
        %figure(f2)
        %subplot(size(c,2),k,(i-1)*k+ord(i,j))
        %showsigcorrs(squeeze(reshape(c{i}(j,:),[r r])),IDX{1},.9)
    end
end


%%
% cd /local_mount/space/dingus/1/RS_analysis/code
% ls
% scp -r param_tunining_analysis.m
% wx2203@s1n1.u19motor.zi.columbia.edu:/home/wx2203/RS_analysis/code
