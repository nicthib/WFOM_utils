%% FC error
clear
load('/local_mount/space/dingus/1/RS_analysis/Draft/params.mat')
addpath(genpath('/local_mount/space/dingus/1/RS_analysis/code'))
DFCopts = [];
DFCopts.nreps = 500;
DFCopts.sig = 0; DFCopts.clustertype = 'kmeans'; 
DFCopts.labelstates = 1; 
DFCopts.donnmf = 0;
DFCopts.dokmeans = 1;
DFCopts.highpassH = 1;
DFCopts.skp = 300;
DFCopts.Hdir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_92comps';
ks = [1 3 5 10 20 30 40];
ws = [1 3 5 10 25 50 100]*20;
runnames = runnames_RS(randperm(numel(runnames_RS),5));
i=1;
%%
for k = 6:numel(ks)
    for w = 6:numel(ws)
        DFCopts.k = ks(k);
        DFCopts.wr = ws(w);
        DFCopts.ww = DFCopts.wr*2+1;
        DFCopts.highpassH = 1;
        [results,~] = getcorrstates_v2(runnames,DFCopts);
        score_n(k,w) = mean(min(results{1}.dj,[],2));
        score_h(k,w) = mean(min(results{1}.dh,[],2));
        DFCopts.highpassH = 0;
        [results,~] = getcorrstates_v2(runnames,DFCopts);
        score_n_nf(k,w) = mean(min(results{1}.dj,[],2));
        score_h_nf(k,w) = mean(min(results{1}.dh,[],2));
        i/(numel(ks)*numel(ws)) 
        i=i+1;
    end
end

%%
close all
yi = interp1(ks,1:numel(ks),7);
xi = interp1(ws/20,1:numel(ws),10);
yinew = interp1(ks,1:numel(ks),7);
xinew = interp1(ws/20,1:numel(ws),25);

figure
subplot(121)
contourf(flipud(rot90(score_n')),0:.05:.5); hold on
scatter(xi,yi,'*','MarkerEdgeColor','r')
yticks(1:numel(ks))
xticks(1:numel(ws))
xticklabels(ws/20)
yticklabels(ks)
caxis([0 .5])
axis image
title('Neural mean correlation distance')
xlabel('Window size (sec)')
ylabel('# of states')

subplot(122)
contourf(flipud(rot90(score_h')),0:.05:.5); hold on
h1 = scatter(xi,yi,'*','MarkerEdgeColor','r');

% Improved model
h2 = scatter(xinew,yinew,'v','MarkerEdgeColor','r');
%
legend([h1 h2],{'Original' 'Proposed'});

yticks(1:numel(ks))
xticks(1:numel(ws))
xticklabels(ws/20)
yticklabels(ks)
yticks(1:numel(ks))
xticks(1:numel(ws))
xticklabels(ws/20)
yticklabels(ks)
caxis([0 .5])
axis image
title('Hemo mean correlation distance')
xlabel('Window size (sec)')
ylabel('# of states')
c = colorbar;
ylabel(c,'Mean corr. distance')
colormap(parula(10))

%save('Hemo_filtered.mat','FC_errh','FC_errn')