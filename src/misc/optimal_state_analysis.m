%% Optimal K for states
% Approach: cross validation
% In this process, the data is partitioned into v parts. Each of the parts 
% is then set aside at turn as a test set, a clustering model computed on 
% the other v − 1 training sets, and the value of the objective function 
% calculated for the test set. 

% This analysis is conditioned off of cm125, cm126 and cm128
clear
%% Behavioral states
root = '/local_mount/space/dingus/1/RS_analysis/H_b_new';
vars = {'H_n','rotf','whisk','pupil'};
% Now get runs with the bahavior vars we want
runnames = getallexps(root,vars)'; 

%% Set options, get runnames for conditioning states
load('TF.mat','cm125')
mouse = 'cm12[568]'; runs = 'BCD';
opts.mapidx = 2; opts.nreps = 50; opts.sig = 6; opts.dokmeans = 1;
opts.ww = 41; opts.skipfactor = (opts.ww-1)/2;
o = cm125.IDXparts{opts.mapidx}; % For visual division of corr maps
divs = [numel(o.L), numel([o.L o.C])]; % For visual division of corr maps
runnames_condition = runnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),runnames)); % Parse out runs that end with B,C,D
% 62 runs, pruned to 60
runnames_condition = runnames_condition(randperm(60,10));

%% Get state centroids
%d=[]; c=[]; cc=[]; st=[];
% Chunk 6, k=12
for k = 20:50
    opts.k=k;
    for ss = 1:10
        tic
        runs_train = runnames_condition;
        idx = ss;%(ss-1)*6+1:ss*6;
        runs_train(idx) = []; % Remove 1/10 test
        runs_test=runnames(idx);
        opts.dokmeans=1; opts.cnt = []; opts.labelstates=0;
        % Get centroid
        [~,c{k,ss},~] = getcorrstates(runs_train,opts);
        % Evaluate distance
        opts.dokmeans=0; opts.cnt=c{k,ss}; opts.labelstates=1;
        [st{k,ss},c{k,ss},d{k,ss},~,cc{k,ss}] = getcorrstates(runs_test,opts);
        disp(['Chunk #' mat2str(ss) ' for k=' mat2str(k) ' done'])
    end
end

%% Get state centroids
%d=[]; c=[]; cc=[]; st=[];
% Chunk 6, k=12
for k = 21:50
    opts.k=k; opts.skipfactor = 10;
    opts.dokmeans=1; opts.labelstates=1;
    [st{k,ss},c{k,ss},d{k,ss},~,cc{k,ss}] = getcorrstates(runnames(1),opts);
    disp(['Chunk #' mat2str(ss) ' for k=' mat2str(k) ' done'])
end

%%
d_ = [];
for k=1:50
    for ss = 1:10
        [d_(k,ss)] = mean(min(d{k,ss}{1},[],2));
    end
end

%%
fitmean = mean(d_(3:end,[1:5,7:10])');
k = 3:50;
r=[]; maxr = 0;
for i = 5:length(fitmean)-5
   l1 = fitmean(1:i);
   l2 = fitmean(i+1:end);
   mdl1 = fitlm(k(1:i),l1);
   mdl2 = fitlm(k(i+1:length(fitmean)),l2);
   r1 = mdl1.Rsquared.Ordinary;
   r2 = mdl2.Rsquared.Ordinary;
   if (r1+r2)/2 > maxr
       maxr = (r1+r2)/2;
       maxmdl1 = mdl1;
       maxmdl2 = mdl2;
   end
end

close all; 
figure
plot(3:50,fitmean)
hold on
plot(maxmdl1)
plot(maxmdl2)
title('Mean label distance (y) vs. k states (x)')
xlabel('# of States')
ylabel('Mean label distance')

figure
plot(3:50,fitmean)
title('Mean label distance (y) vs. k states (x)')
xlabel('# of States')
ylabel('Mean label distance')
