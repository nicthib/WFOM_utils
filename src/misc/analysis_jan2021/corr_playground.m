%% Initialize
clear
Hdir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; % The base folder we are pulling data from
vars = {'H'}; % Variables required to include runs in our analysis
allrunnames = getallexps(Hdir,vars)'; % Gets all runs in the folder with the criteria above

%% Set options
%clearvars -EXCEPT allrunnames
opts.k = 4; opts.nreps = 100; % Some options
opts.ww = 41; opts.sig = 6; opts.dokmeans = 1;
opts.skipfactor = (opts.ww-1)/2; % Some options
opts.Hdir = Hdir;
runs = 'BCD'; mouse = 'cm12[5]_[3]';
runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames)); % Parse out runs that end with B,C,D
runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames)); % Parse out runs that are 'mouse'

%% Get states
opts.variable = 'jrgeco';
[st,c,d,~,~,p] = getcorrstates(runnames,opts);

%% Get onsets
runningBlocks = [];
for i = 1:numel(runnames)
   load(fullfile(Hdir,runnames{i})) % Load data 
   runningBlocks = [runningBlocks getrunningblocks(H.jrgeco,m,400,d{i},p{i})];
end

%% Visualize states w/ behavior
n = 1; m = 1;
svb_vis(runnames{n},st{m}{n},c{m},divs,opts.mapidx,[]);

%% Visualize state centroids
show_state_centroids(c{1},[9 16 9 5 7],1:5);

%% Transition matrix
%subplot(2,k,k+1:k+3)
close all
t_mat = gettransitionmatrix(cat(1,st{1}{:}),1:k,1);
G = [];
for i = 1:size(t_mat,1)
    for j = 1:size(t_mat,1)
        if t_mat(j,i) > 0
            G(end+1,:) = [i j t_mat(i,j)];
        end
    end
end
g = digraph(G(:,1),G(:,2),G(:,3));
L = 5*g.Edges.Weight/max(g.Edges.Weight);
%subplot(2,k,k+4:k+7)
plot(g,'LineWidth',L);
axis off
hold on

%%

for i = 1:1000
    st1 = st{1}{1}(i);
    idx = min(find(st{1}{1}(i:end) ~= st1)); idx = idx+i-1;
    st2 = st{1}{1}(idx);
    rat = d{1}{1}(i,st1)/d{1}{1}(i,st2);
    xc(i) = (x(st1)*rat + x(st2)*(1-rat));
    yc(i) = (y(st1)*rat + y(st2)*(1-rat));
end

%% Video
vo.frames = [0 4000];
vo.timestamp = 1;
vo.mouse = 'cm125_6';
vo.run = 'runC';
[vid,~] = makeFLIRvid(vo);

%% Other video
load(runnames{n}{m})
slcorr_video(H_n{opts.mapidx},st{n}{m},c{n},getwin(opts.ww,opts.sig),vid,strrep(vid,'_.avi','_tile.avi'));

%%
[~,x,y] = makeFLIRvid(mouse,'B','1',[idxs(1)-30 idxs(1) + 90]);
opts.x = x; opts.y = y; opts.crop = 0;
for i = 2:36
    makeFLIRvid(mouse,'B',mat2str(i),[idxs(i)-30 idxs(i)+90],opts);
end

%% making cool video?
% State 1 is of interest
idx = strfind(st1',[1 5]);
trans_15 = [];
for i = 1:numel(idx)
    if max(rot(idx(i)-30:idx(i)+30)) < 2
        trans_15 = [trans_15 idx(i)];
    end
end
idxs = trans_15(1:36);

%% Body movement score
%clear
runnames = 'BCDEFGHIJ';
sessionname = 'cm125_3';
B_TC_full = [];
p_vec_full = [];
rotf_full = [];
for i = 1:numel(runnames)
    tmp = extract_body_movement_score([sessionname '_run' runnames(i)]);
    load([sessionname '_run' runnames(i) '_H.mat'],'m','rotf')
    try rotf = m.rotf; catch; end
    rotf_full = [rotf_full rotf];
    B_TC_full = [B_TC_full tmp(1:11980)];
end
% Now do states above

%%
close all;
d_full = cell2mat(d');
p_vec_full = cell2mat(p_vec');
tiledlayout(3,1);
ax(1) = nexttile;
plot(B_TC_full); hold on
plot(abs(rotf_full/max(rotf_full)))
legend('Body movement score','Locomotion')

ax(2) = nexttile;
plot(d_full(:,1)); hold on
plot(d_full(:,12))
legend('State 1 distance','State 12 distance')


ax(3) = nexttile;
plot(p_vec_full(:,1)); hold on
plot(p_vec_full(:,12))
legend('State 1 distance','State 12 distance')


linkaxes(ax,'x')
pan xon