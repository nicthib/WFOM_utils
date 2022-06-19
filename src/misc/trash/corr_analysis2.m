%% Correlation analysis with dist as a function of span
% Setup, options
clearvars; %close all
addpath /local_mount/space/revault/revault2/Ewoud/SRGAP2C/Colormap/
addpath(genpath('/local_mount/space/dingus/1/RS_analysis/'))
root = '/local_mount/space/dingus/1/RS_analysis/H_b';
load ERES_2.mat
c = [];
opts.mapidx = 8; opts.k = 7; opts.nreps = 30; opts.fsp = 10; 
opts.skipfactor = 10; opts.ww = 41; opts.sig = 0; opts.dokmeans = 1;

%%
% %expname = {'cm125_5','cm125_6','cm126_5','cm126_6','cm128_4','cm128_6'};
% expname = {'cm125_5'};
% runs = 'BCD'; %HIJ';
% runnames = {};
% for i = 1:numel(expname)
%     for j = 1:numel(runs)
%         runnames{end+1} = [expname{i} '_run' runs(j) '_H.mat'];
%     end
% end

exp_root = '/local_mount/space/dingus/1/RS_analysis/H_b';
expname = {'cm125_5','cm125_6'};%,'cm126_5','cm126_6','cm128_4','cm128_6'};
expname = {'cm128_3'};
runs = 'D';
runnames = {};

for i = 1:numel(expname)
    for j = 1:numel(runs)
        runnames{end+1} = fullfile(exp_root,[expname{i} '_run' runs(j) '_H.mat']);
    end
end


%%

[st,c,d] = getcorrstates(runnames,opts);


%% Visualize
load(fullfile(root,'TF.mat'),'cm125')
o = cm125.IDXparts{opts.mapidx}; % For visual division of corr maps
divs = [numel(o.L), numel([o.L o.C])]; % For visual division of corr maps
n = 6;
svb_vis(runnames{n},st{n},c,divs,opts.mapidx);




%%
% n = 1;
% % for a = 1:numel(runnames)
%     for b = 1:numel(winsizes)
%         opts.sp = winsizes(b)*10;
%         %[st{n},c{n},cc{n},IDX{n},spk{n},d{n}] = getcorrstates(runnames{a},opts);
%         [st{n},c{n},d{n}] = getcorrstates(runnames,opts);
%         n = n+1;
%     end
% % end

[st,c,d] = getcorrstates(runnames,opts);



% 
% r = sqrt(size(c,2));
% k = opts.k;
% 
% ord = repmat(1:k,[size(c,2) 1]);
% for i = 2:size(c,2)
%     c_dist = pdist2(c{1},c{i},'euclidean');
%     ordtmp = zeros(1,k);
%     while sum(ordtmp > 0) < k
%         [nmx,nmy] = find(c_dist==min(c_dist(:)));
%         if ordtmp(nmy) == 0 && isempty(find(ordtmp==nmx))
%             ordtmp(nmy) = nmx;
%         end
%         c_dist(nmx,nmy) = NaN;
%     end
%     ord(i,:) = ordtmp;
% end
% st2 = {};
% for j = 2:numel(st)
%     for i = 1:k
%         st2{j}(st{j}==i) = ord(j,i);
%     end
% end
% st2{1} = st{1}';
% disp('DONE!')

%%
% vO = VideoWriter('Corr_state_comparison.avi');
% vO.FrameRate = 20;
% open(vO)
close all
figure
dists = [];
for t = 1:1000
    subplot(3,3,2)
    for i = 1:numel(winsizes)
        subplot(3,3,i)
        imagesc(reshape(cc{i}(t,:),[sqrt(676) sqrt(676)]))
        caxis([0 1])
        title([mat2str(winsizes(i)) ' sec window'])
        axis image;  colormap jet; set(gca,'Xticklabel',[]); set(gca,'Yticklabel',[])
        subplot(3,3,i+3)
        imagesc(reshape(c{i}(st{i}(t),:),[sqrt(676) sqrt(676)]))
        caxis([0 1])
        axis image;  colormap jet; set(gca,'Xticklabel',[]); set(gca,'Yticklabel',[])
        subplot(3,3,i+6)
        dists(i,t) = pdist2(cc{i}(t,:),c{i}(st{i}(t),:));
        plot((1:t)/20,dists(i,:)); hold off
        xlim([t-100 t]/20)
        ylim([0 20])
        subplot(331)
        ylabel('dFC map')
        subplot(334)
        ylabel('State nearest')
        subplot(337)
        ylabel('label/corr distance')
        subplot(338)
        xlabel('time (sec)')
    end
%     writeVideo(vO,getframe(gcf))
end
% close(vO)


%% State consistency bar plot
clear st_*
st2mat = cell2mat(st2')';

for i = 1:size(st2mat,2)
    for j = 1:k
        st_f(i,j) = numel(find(st2mat(:,i)==j))/size(cc{1},1);
    end
end

% Table for reprated ANOVA
inds = [1:3];
inds_pre = [inds inds+6 inds+12 inds+18 inds+24 inds+30];
inds = [4:6];
inds_post = [inds inds+6 inds+12 inds+18 inds+24 inds+30];
labels = cell(1,36);
labels(:) = {'pre'};
labels(inds) = {'post'};
t = table(labels',st_f(:,1),st_f(:,2),st_f(:,3),st_f(:,4),st_f(:,5),st_f(:,6),st_f(:,7)...
    ,'VariableNames',{'condition','k1','k2','k3','k4','k5','k6','k7'});
%%
close all
st_fm(:,1) = mean(st_f(inds_pre,:),1);
st_fm(:,2) = mean(st_f(inds_post,:),1);

st_fs(:,1) = std(st_f(inds_pre,:),[],1);
st_fs(:,2) = std(st_f(inds_post,:),[],1);

bar(1:k,st_fm)
hold on
er = errorbar(1:k,st_fm,st_fs,st_fs);
er.Color = [0 0 0];
er.LineStyle = 'none';







