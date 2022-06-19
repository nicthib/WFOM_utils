% [st,cnt,d,dmin,cc,p] = getcorrstates(runnames,opts) calculates 
% correlation states and associated variables for a list of given datasets.
%
% Inputs:
% runnames: cell array of runname scripts. These correspond to datasets
% that contain the cell matrix H_n, with timecourses that correspond to
% neural timecourses.
% 
% opts: structure containing options for calculating correlation states:
%
% opts.k: Number of states to calculate.
%
% opts.ww: The correlation window width in # of frames. This number must be
% odd.
%
% opts.sig: setting this to a value other than 0 tapers the window by a
% gaussian window with the given sigma.
%
% opts.nreps: # of replicates used for k-means algorithm.
% 
% opts.skipfactor: the the skip value for the sliding window. This saves
% computational time during k-means. A generally good value is abotu half
% the window width.
%
% opts.dokmeans: if set to 0, k-means will not be performed, and outputs
% will be calculated using previously calculated states (opts.cnt)
%
% opts.variable: set to whatever fieldname of H you want to use for state
% analysis
%
% Outputs: 
% st: a vector containing the state label of every frame from the input
% dataset. Note that opts.skipfactor is ignored here: st will contain state
% values for the entire dataset.
%
% cnt: state centroids. This will take the shape k x n^2, where k = opts.k,
% and n is the number of rows in the H_n dataset.
%
% d: similar to st, but contains the raw distances for each state from the
% raw correlation matrix. Note that argmin(d(i)) = st(i)
%
% cc: full corr matrix

function [results, cnt] = getcorrstates_v2(runnames,opts)

% Set default options
if ~isfield(opts,'dokmeans'); opts.dokmeans = 1; end
if ~isfield(opts,'labelstates'); opts.labelstates = 1; end
if ~isfield(opts,'variable'); opts.variable = 'jrgeco'; end
if ~isfield(opts,'nreps'); opts.nreps = 1; end
if ~isfield(opts,'Hdir'); opts.Hdir = []; end
if ~isfield(opts,'calc_prob'); opts.calc_prob = 1; end

results = [];
cc_ds = [];
cnt = [];

warning('off','all') % K-means tends to send warnings while in parallel.

% Compose correlation matrices
opts.win = getwin(opts.ww,opts.sig); % Create window (square if sig = 0, tapered otherwise)
h = waitbar(0,['Getting sliding correlations for ' num2str(numel(runnames)) ' runs...']); h.Position(1:2) = [0 1080];
for i = 1:numel(runnames) % iterate through each cell of H.
    load(fullfile(opts.Hdir,runnames{i}),'H','m') % Load data
    Htmp = H.(opts.variable);
    subregions = [9 16 9 5 7];
    subregions = [0 cumsum([subregions subregions])];
    Htmp2 = [];
    for j = 1:numel(subregions)-1
       Htmp2(j,:) = mean(Htmp(subregions(j)+1:subregions(j+1),:),1);
    end
    Htmp = Htmp2;
    slc_table = struct2table(slcorr_table(Htmp,m,opts));
    cc_ds{i} = slc_table.cc;    
    if opts.labelstates % If labelling states, get the full corr matrix
        opts_tmp = opts; opts_tmp.skp = 1;
        slc_table = struct2table(slcorr_table(Htmp,m,opts_tmp));
        results{i}.cc = slc_table.cc;
    end
    waitbar(i/numel(runnames),h)
end
close(h)
cc_ds = cell2mat(cc_ds');
%

% Perform clustering
options = statset('UseParallel',1); % To increase speed
if opts.dokmeans == 1
    disp(['Performing k-means on ' mat2str(size(cc_ds,1)) ' windows...'])
    if ~isfield(opts,'cnt')
        [st_tmp] = kmeans(cc_ds,opts.k,'Distance','correlation','options',options,'Replicates',opts.nreps);
    else
        [st_tmp] = kmeans(cc_ds,opts.k,'Distance','correlation','options',options,'Start',opts.cnt);
    end
    
    for i = 1:opts.k
        cnt(i,:) = mean(cc_ds(find(st_tmp==i),:),1);
    end
else
    cnt = opts.cnt;
end
%

% Label states
if opts.labelstates
    disp(['Labeling states for ' num2str(numel(results)) ' runs...']);
    for i = 1:numel(results)
        results{i}.d = pdist2(results{i}.cc,cnt,'correlation'); % label each full matrix using output centroids
        [~,results{i}.st] = min(results{i}.d,[],2); % state is determined by minimum dist to centroids.
    end
end

disp('Done')

% Determine state probabilities
% if opts.calc_prob
%     disp(['Getting state probabilities for ' num2str(numel(results)) ' runs...']);
%     p = [];
%     for i = 1:numel(results)
%         for j = 1:size(results{i}.cc,1)
%             V1 = results{i}.cc(j,:);
%             V1 = V1/norm(V1);
%             for k = 1:opts.k
%                 V2 = cnt(k,:)';
%                 V2 = V2/norm(V2);
%                 results{i}.p(j,k) = (V1 * V2)^2;
%             end
%         end
%     end
% end
