%% How do states compare (neural vs. hemo)?
%% Initialize
clear
H_dir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; % The base folder we are pulling data from
vars = {'H'}; % Variables required to include runs in our analysis
allrunnames = getallexps(H_dir,vars)'; % Gets all runs in the folder with the criteria above

%% Set options
%clearvars -EXCEPT allrunnames
opts.k = 7; opts.nreps = 100; % Some options
opts.ww = 201; opts.sig = 23; opts.dokmeans = 1;
opts.skipfactor = (opts.ww-1)/2; % Some options
runs = 'BCD'; mouse = 'cm12[5678]_[3456]';
runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames)); % Parse out runs that end with B,C,D
runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames)); % Parse out runs that are 'mouse'

%% Get states
opts.variable = 'jrgeco';
[st{1},c{1},d{1},~,~,p{1}] = getcorrstates(runnames,opts);
opts.variable = 'chbr';
[st{2},c{2},d{2},~,~,p{2}] = getcorrstates(runnames,opts);

%% Brute force st matching
pms = perms(1:opts.k);
d1 = d{1}{1};
s = zeros(size(pms,1),1);
for i = 1:size(pms,1)
   d2 = d{2}{1}(:,pms(i,:));
   s(i) = sum(sum(abs(d2-d1)));
end
bestord = pms(find(s==min(s)),:);

%%
close all
n=2;
tiledlayout(3,1)
load(fullfile(H_dir,runnames{n}),'m')
ax(1) = nexttile;
imagesc([d{1}{n}]')
caxis([0 1])
xlabel('frame #'); ylabel('State #')
title('Neural state distance')
ax(2) = nexttile;
imagesc([d{2}{n}(:,bestord)]')
caxis([0 1])
xlabel('frame #'); ylabel('State #')
title('chbr state distance')
ax(3) = nexttile;
wr = (opts.ww-1)/2;
plot([m.whisk(wr:end) zeros(1,wr)]); hold on; plot([m.rotf(wr:end) zeros(1,wr)])
linkaxes(ax,'x')
colormap gray
legend('whisking','locomotion')

%%
close all
show_state_centroids(c{1},[9 16 9 5 7],1:5);
show_state_centroids(c{2}(bestord,:),[9 16 9 5 7],1:5);
