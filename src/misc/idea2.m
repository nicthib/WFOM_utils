%% Initialize
clear
Hdir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; % The base folder we are pulling data from
fundir = '/local_mount/space/dingus/1/RS_analysis/code';
addpath(genpath(fundir));
load('/local_mount/space/dingus/1/RS_analysis/analysis_feb2021/allrunnames.mat')
%%
opts.k = 7; opts.nreps = 100; % Some options
opts.ww = 31; opts.sig = 6; opts.dokmeans = 1;
opts.skipfactor = 100;%(opts.ww-1)/2; % Some options
opts.Hdir = Hdir;
opts.variable = 'jrgeco';
runs = 'BCD'; mouse = 'cm12[5]_[4]';
runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames)); % Parse out runs that end with B,C,D
runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames)); % Parse out runs that are 'mouse'

%% Get states
[st{1},c{1},d{1},~,~,p{1}] = getcorrstates(runnames,opts);

%% Get bootstrapped TMs, eigenvectors for stability analysis
st_all = cell2mat(st{1}');
st1 = st_all(1:end-1); st2 = st_all(2:end);
Vd = [];
for i = 1:1000
    idx = randperm(numel(st_all)-1,round(numel(st_all)/10));
    T = gettransitionmatrix(st1(idx),st2(idx),opts.k,0);
    [V,D] = eig(T/sum(T(:))); Vd(:,i) = diag(V);
    Td(:,i) = diag(T)/20;
end

%% Show
close all
for i = 1:opts.k
    errorbar(i,mean(Td(i,:)),std(Td(i,:))*2)
    hold on
end
xlim([.8 opts.k+.2])
title('Stability values for neural states (eigenvector of 100 bootstrapped TMs')
ylabel('Stability value')
xlabel('State #')
%ylim([0 .2])

%% Transition matrix
T = gettransitionmatrix(st1,st2,opts.k,0);
imagesc(log10(T/sum(T(:))))
axis image
colormap hot
colorbar
caxis([-5 -2])
xlabel('state at t+1')
ylabel('state at t')
title('State transition matrix (log scale)')

%% Direct measurement of state lengths

for i = 1:opts.k
    tmp = st_all==i;
    tmp(1) = 0;
    tmp(end) = 0;
    tmpd = diff(tmp);
    st_start = find(tmpd == 1);
    st_end = find(tmpd == -1);
    stl = st_end-st_start;
    errorbar(i,mean(stl),std(stl))
    hold on
end