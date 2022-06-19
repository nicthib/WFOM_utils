
%% START HERE
clear; close all
Hdir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; % The base folder we are pulling data from
fundir = '/local_mount/space/dingus/1/RS_analysis/code';
addpath(genpath(fundir));
load('/local_mount/space/dingus/1/RS_analysis/Draft/params.mat','runnames_singlemouse_goodFC','final_ord','m_all')
savedir = '/local_mount/space/dingus/1/RS_analysis/Draft/part5';
IDX = remap_IDX(m_all.cm125.IDX_final,final_ord); clear m_all

%% Set options
opts = [];
opts.k = 7; opts.nreps = 500; % Some options
opts.ww = 201; opts.sig = 0; 
opts.dokmeans = 1; opts.labelstates = 1;
opts.skp = 1000; % Some options
% for i = 1:5
%     runnames_singlemouse{i} = runnames_RS(cellfun(@(s) ~isempty(regexp(s,['cm12' mat2str(i+3)])),runnames_RS)); % Parse out runs that are 'mouse'
% end
runnames_singlemouse = runnames_singlemouse_goodFC;
runnames_singlemouse{4}(1:3) = [];
runnames_label = [runnames_singlemouse{1};...
                  runnames_singlemouse{2};...
                  runnames_singlemouse{3};...
                  runnames_singlemouse{4}
                  runnames_singlemouse{5}];

% cm127 fix?

%% Get states
opts.labelstates = 0;
[~,cnt_label] = getcorrstates_v2(runnames_label,opts);   
opts.labelstates = 1;
[~,I] = sort(mean(cnt_label,2)); opts.cnt = cnt_label(I,:);
opts.dokmeans = 0;
for n = 1:numel(runnames_singlemouse)
    [results{n},cnt{n}] = getcorrstates_v2(runnames_singlemouse{n},opts);
end

%% Show dFC process
%[results_tmp,~] = getcorrstates_v2(runnames_singlemouse{n}(1),opts);
close all
load('/local_mount/space/dingus/1/RS_analysis/Draft/params.mat', 'final_ord')
networks = [8 8 7 4 7 12];
t = [300 400 500];
load(runnames_singlemouse{n}{1},'H')
figure('Position',[675 625 1246 336])
imagesc(H.jrgeco(final_ord,1:1000)); 
caxis([-.05 .05]); colormap coolwarm; colorbar; axis image; axis off


figure
for i = 1:3
    subplot(1,3,i)
    cc_tmp = reshape(results_tmp{1}.cc(t(i),:),[92 92]);
    cc_tmp = cc_tmp(final_ord,final_ord);
    imagesc(cc_tmp); 
    hold on; axis image; axis off;
    add_network_lines(networks)
    caxis([0 1])
    title(['Win center = ' mat2str(t(i)/20)])
end
colormap jet
hgexport(1,fullfile(savedir,'5A.eps'))
hgexport(2,fullfile(savedir,'5B.eps'))


%% Load behavior
load('/local_mount/space/dingus/1/RS_analysis/Draft/part3/behavior.mat')
for n = 1:numel(results)
    for i = 1:numel(results{n})
        fn = fieldnames(behavior{n}{i});
        for f = 1:numel(fn)
            results{n}{i}.(fn{f}) = behavior{n}{i}.(fn{f});
        end
    end
end

%% Get onsets (only need to run once)
T1 = [300 500];
for n = 1:numel(results)
    runningBlocks{n} = [];
    for i = 1:numel(results{n})
        try runningBlocks{n} = [runningBlocks{n} getrunningblocks(results{n}{i},T1,runnames_singlemouse{n}{i})]; 
        catch; disp('failed'); end
    end
end
disp('Done')

%% prune
n = 1;
Tr = 200; % # frames rest period (higher = less N)
runningBlocksTable = [];
for i = 1:4
    runningBlocksTable = [runningBlocksTable; sortrows(struct2table(runningBlocks{i}), 'offsetduration', 'descend')];
end

TableForAnalysis = runningBlocksTable;
TableForAnalysis(TableForAnalysis.offsetduration < Tr,:) = [];
TableForAnalysis(TableForAnalysis.onsetduration < 50,:) = [];
TableForAnalysis(TableForAnalysis.p_score < .9,:) = [];

%%
tp = .2;
bs = [];
close all
mice = {'cm125','cm126','cm127','cm128'};
d_all = []; p_all = []; r_all = []; w_all = [];
for i = 1:size(TableForAnalysis,1)
    l = TableForAnalysis.onset(i);
    idx = 1:sum(T1);
    d_all(:,:,end+1) = TableForAnalysis.d{i}(idx,:);
    p_all(:,end+1) = TableForAnalysis.pupil{i}(idx,:);
    r_all(:,end+1) = TableForAnalysis.rotf{i}(idx,:);
    w_all(:,end+1) = TableForAnalysis.whisk{i}(idx,:);
end
%d_all = d_all*100;

colormap jet
cmap = lines(opts.k);
t = linspace(-T1(1)/20,T1(2)/20,size(d_all,1));
for i = 1:opts.k
    subplot(211)
    bl = squeeze(mean(d_all(1,i,:)));
    stdshade(squeeze(d_all(:,i,:))'-bl,tp,cmap(i,:),t+opts.ww/(20*2));
    hold on
end
xlim([-10 25])
subplot(212)
hold on
B_cmap = lines(3);
B_sf = [1 .5 .2];
stdshade((w_all'-nanmean(w_all(1,:),2))/B_sf(1),tp,B_cmap(1,:),t);
stdshade((r_all'-nanmean(r_all(1,:),2))/B_sf(2),tp,B_cmap(2,:),t);
stdshade((p_all'-nanmean(p_all(1,:),2))/B_sf(3),tp,B_cmap(3,:),t);
%ylim([-.2 1.5])    
xlim([-10 25])

% delete shades
subplot(211)
axesHandlesToChildObjects = findobj(gca, 'Type', 'line');
if ~isempty(axesHandlesToChildObjects)
    delete(axesHandlesToChildObjects);
end  
axis off
subplot(212)
axesHandlesToChildObjects = findobj(gca, 'Type', 'patch');
if ~isempty(axesHandlesToChildObjects)
    delete(axesHandlesToChildObjects);
end  
hgexport(1,fullfile(savedir,'6A_shades.eps'))
axis off

%% Sig corrs
close all
cnt_show = reshape(opts.cnt,opts.k,92,92);
cnt_show = cnt_show(:,final_ord,final_ord);
IDX = remap_IDX(m.IDX,final_ord);
for i = 1:7
    figure(1)
    subplot(1,7,i)
    imagesc(squeeze(cnt_show(i,:,:))); hold on
    add_network_lines(networks)
    axis image; axis off; caxis([0 1]); colormap jet

    figure(2)
    subplot(1,7,i)
    showsigcorrs(squeeze(cnt_show(i,:,:)),IDX,.8,8000)
end

hgexport(1,fullfile(savedir,'6B-1.eps'))
hgexport(2,fullfile(savedir,'6B-2.eps'))

%% Example run
t = linspace(0,600,size(results{n}{m}.d,1));
close all; 
figure
n=1; m=10;
plot(t,results{n}{m}.whisk); hold on
plot(t,results{n}{m}.rotf);
plot(t,results{n}{m}.pupil-15);
legend('Whisk','Movement','Pupil','Location','EastOutside')
xlim([1000 3500]/20)

figure
for i = 1:7
    plot(t,results{n}{m}.d(:,i));% + repmat([0;1;2;3;4;5;6]/5,[1 size(results{n}{m}.d,1)])')
    hold on
end
legend('State 1','State 2','State 3','State 4','State 5','State 6','State 7','Location','EastOutside')    
xlim([1000 3500]/20-5)
xlabel('Time (sec)')


figure
imagesc(results{n}{m}.st')
xlim([1000 3500]-100)
colormap(lines(7))
caxis([0 8])
axis off

hgexport(1,fullfile(savedir,'6C.eps'))
hgexport(2,fullfile(savedir,'6D.eps'))
hgexport(3,fullfile(savedir,'6E.eps'))

