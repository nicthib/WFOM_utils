%% High level analysis for RS group
%% 6/18/2020
% This data is conditioned off of cm125, cm126 and cm128
clear
%% Behavioral states
root = '/local_mount/space/dingus/1/RS_analysis/H_b_new';
vars1 = {'H_n','rotf','whisk','pupil'};
vars2 = {'H_n'};
% First get as many good runs as possible
runnames_all = getallexps(root,vars2)'; 
% Now get runs with the bahavior vars we want
runnames_behav = getallexps(root,vars1)'; 

%% Set options, get runnames for conditioning states
load('TF.mat','cm125')
mouse = 'cm12[8]'; runs = 'BCD';
opts.mapidx = 3; opts.k = 20; opts.nreps = 50; opts.sig = 6; opts.dokmeans = 1; opts.dobehavior = 0;
opts.ww = 41; opts.skipfactor = (opts.ww-1)/2;
o = cm125.IDXparts{opts.mapidx}; % For visual division of corr maps
divs = [numel(o.L), numel([o.L o.C])]; % For visual division of corr maps
runnames_condition = runnames_all(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),runnames_all)); % Parse out runs that end with B,C,D
runnames_condition = runnames_condition(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames_condition)); % Parse out runs that are 'mouse'

%% Get state centroids
opts.ww = 41; opts.skipfactor = (opts.ww-1);
[~,c.w2,~] = getcorrstates(runnames_condition,opts);
%opts.ww = 101; opts.skipfactor = (opts.ww-1);
%[~,c.w5,~] = getcorrstates(runnames_condition,opts);
%opts.ww = 401; opts.skipfactor = (opts.ww-1);
%[~,c.w20,~] = getcorrstates(runnames_condition,opts);

%% Get state variables
opts.dokmeans = 0;
opts.cnt = c.w2; opts.ww = 41; [st.w2,~,~] = getcorrstates(runnames_behav,opts);
opts.cnt = c.w5; opts.ww = 101; [st.w5,~,~] = getcorrstates(runnames_behav,opts);
opts.cnt = c.w20; opts.ww = 401; [st.w20,~,~] = getcorrstates(runnames_behav,opts);

%% Get behavioral variables 
stb = getbehaviorz(runnames_behav);

%% Labelling one mouse using general states


