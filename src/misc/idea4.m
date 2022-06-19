% idea 4: look at state behavior as a funtion of time since/till next
% movement

%% Initialize
clear
Hdir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; % The base folder we are pulling data from
fundir = '/local_mount/space/dingus/1/RS_analysis/code';
addpath(genpath(fundir));
load('/local_mount/space/dingus/1/RS_analysis/analysis_feb2021/allrunnames.mat')

%% Set options
opts.k = 7; opts.nreps = 100; % Some options
opts.ww = 61; opts.sig = 0; opts.dokmeans = 1;
opts.skipfactor = (opts.ww-1)/2; % Some options
opts.Hdir = Hdir; opts.variable = 'jrgeco';
runs = 'BCD'; mouse = 'cm12[5]_[4]';
runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames)); % Parse out runs that end with B,C,D
runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames)); % Parse out runs that are 'mouse'

%% Get states
[st,c,d,~,~,p] = getcorrstates(runnames,opts);

%% Get behavioral info
for i = 1:numel(runnames)
    load(fullfile(Hdir, runnames{i}),'m')
    rz = getbehavioralbouts(m.rotf,.5,10,10);
    wz = getbehavioralbouts(m.whisk,.5,10,10);
    [ftn_r{i},ftp_r{i}] = behav_time(rz);
    [ftn_w{i},ftp_w{i}] = behav_time(wz);
end
d_all_vec = cell2mat(d_all);
st_all_vec = cell2mat(st_all);

%% Look at time since last movement
ftp_r_all = cell2mat(ftp_w');
ftn_r_all = cell2mat(ftn_w');
for i = 1:500
    %idx = find(and(ftp_r_all==i,ftn_r_all > 500));
    idx = find(ftp_r_all==i);
    d_mean(:,i) = mean(d_all_vec(idx,:),1);
    for j = 1:opts.k
        st_mean(j,i) = sum(st_all_vec(idx)==j)/numel(idx);
    end
    n(i) = numel(idx);
end

%%
close all
cmap = lines(opts.k);
plotvec = 1:min(find(n<=100));
for i = 1:opts.k
   %plot(smooth(d_mean(i,:),1000),'Color',cmap(i,:))
   hold on
   plot(linspace(0,plotvec(end)/20,plotvec(end)),d_mean(i,plotvec),'Color',cmap(i,:))
   legend
end
xlabel('time (sec)')