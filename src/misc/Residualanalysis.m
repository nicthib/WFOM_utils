%% Initialize
clear
root = '/local_mount/space/dingus/1/RS_analysis/H_b';
vars = {'H_n'};
allrunnames = getallexps(root,vars)'; % Gets all runs in the H_nf folder
load('TF.mat','cm125')
opts.mapidx = 3; opts.k = 7; opts.nreps = 500; opts.sig = 6;
opts.ww = 41; opts.skipfactor = (opts.ww-1)/2;
%% Get states
mice = {'cm','cm125','cm126','cm127','cm128'}; runs = 'BCD';
for i = 1:numel(mice)
opts.dokmeans = 1;
runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames))'; % Parse out runs that end with B,C,D
runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mice{i})),runnames))'; % Parse out runs that are 'mouse'

[~,c1,~] = getcorrstates(runnames,opts);
mouse = 'cm'; runs = 'BCD';
runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames))'; % Parse out runs that end with B,C,D
runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames))'; % Parse out runs that are 'mouse'
opts.dokmeans = 0; opts.cnt = c1;
[st{i},c{i},d{i}] = getcorrstates(runnames,opts);
end

%%

rndivs = {1:17,18:31,32:40,41:60};
names = {'cm125','cm126','cm127','cm128'};
d_a = [];
close all; figure; hold on
cmap = [0 0 0;1 0 0;0 0 1];
for i = 1:numel(d{1})
    d_a(i) = mean(min(d{1}{i},[],2));
end
for i = 1:4
    dtmp = [];
    for j = 1:numel(d{i+1})
        dtmp(j) = mean(min(d{i+1}{j},[],2));
    end
    
    
    errorbar(i-.2,mean(d_a(rndivs{i})),std(d_a(rndivs{i})),'Color',[0 0 0])
    errorbar(i-.1,mean(dtmp(rndivs{i})),std(dtmp(rndivs{i})),'Color',[1 0 0])
    errorbar(i-.0,mean(dtmp(tmp)),std(dtmp(tmp)),'Color',[0 0 1])
    tmp = 1:4; tmp(i) = [];
    tmp = cell2mat(rndivs(tmp));
    scatter(repmat(i-.2,[numel(rndivs{i}) 1])',d_a(rndivs{i}),'MarkerEdgeColor',cmap(1,:))
    scatter(repmat(i-.1,[numel(rndivs{i}) 1])',dtmp(rndivs{i}),'MarkerEdgeColor',cmap(2,:))
    scatter(repmat(i-.0,[numel(tmp) 1])',dtmp(tmp),'MarkerEdgeColor',cmap(3,:))

end
legend('States obtained from all mice','States obtained from this mouse','States obtained from another mouse')
xlim([0 5])
ylim([0 1000])
set(gca,'XTick',[1:4],'xticklabel',names)
ylabel('Average label distance from centroid')
title('Label distance statistics (states obtained from cm128)')

%%

rndivs = {1:17,18:31,32:40,41:60};
names = {'cm125','cm126','cm127','cm128'};

d_a = [];
close all; figure; hold on
cmap = [0 0 0;1 0 0;0 0 1];
n = 1;
for i = 1:numel(d{n})
    d_a(i,1) = mean(min(d{n}{i},[],2));
    mlist(i,:) = runnames{i}(1:5);
end
violinplot(d_a,cellstr(mlist));
xlim([0 5])


