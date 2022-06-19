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
% dmin: column-wise minimum from d.
%
% cc: full corr matrix
%
% p: state probability, calculated as innter cosine product.

function [st,cnt,d,dmin,cc,p] = getcorrstates(runnames,opts)

% Set default options
if ~isfield(opts,'dokmeans'); opts.dokmeans = 1; end
if ~isfield(opts,'labelstates'); opts.labelstates = 1; end
if ~isfield(opts,'calc_prob'); opts.calc_prob = 1; end
if ~isfield(opts,'variable'); opts.variable = 'jrgeco'; end
if ~isfield(opts,'nreps'); opts.nreps = 1; end
if ~isfield(opts,'Hdir'); opts.Hdir = []; end
if ~isfield(opts,'running'); opts.running = 'both'; end
%

% Compose H matrices
warning('off','all') % K-means tends to send warnings while in parallel.
Hall = [];  d = []; st = []; % Initialize
for i = 1:numel(runnames)
    load(fullfile(opts.Hdir,[runnames{i} '_H.mat']),'H','m') % Load data
    Htmp = H.(opts.variable); % Make sure each H is the same # of frames
    Htmp(isnan(Htmp)) = 0; Htmp(isinf(Htmp)) = 0;
    Hall{end+1} = Htmp;
end
%

% Compose correlation matrices
win = getwin(opts.ww,opts.sig); % Create window (square if sig = 0, tapered otherwise)
cc_ds = []; % cc_ds is the temporally downsampled corr matrix for kmeans (to save computation time).
cc = []; % cc is the full corr matrix for labeling.
h = waitbar(0,['Getting sliding correlations for ' num2str(numel(Hall)) ' runs...']); h.Position(1:2) = [0 1080];
for i = 1:numel(Hall) % iterate through each cell of H.
    if opts.dokmeans % If doing kmeans, get the downsampled corr matrix
        cc_ds{i} = atanh(slcorr(Hall{i},opts.skp,win));%,m.rotf)); % ADDED ATANH() 2/5/21
    end
    if opts.labelstates % If labelling states, get the full corr matrix
        cc{i} = atanh(slcorr(Hall{i},1,win));
    end
    waitbar(i/numel(Hall),h)
end
close(h)
cc_ds = cell2mat(cc_ds)'; % Concatenate cells, then permute to cluster temporally and not spatially
%

% Perform clustering
if opts.dokmeans
    options = statset('UseParallel',0); % To increase speed
    disp('Performing k-means...')
    if isfield(opts,'cnt')
        [st_tmp,cnt] = kmeans(cc_ds,opts.k,'Distance','correlation','options',options,'Start',opts.cnt);
    else
        [st_tmp,cnt] = kmeans(cc_ds,opts.k,'Distance','correlation','Replicates',opts.nreps,'options',options);
    end
    % Recapitulate centroids
    for i = 1:opts.k
        cnt(i,:) = tanh(mean(cc_ds(find(st_tmp==i),:),1));
    end
    [~,I] = sort(sum(cnt,2)); % sort by centroid sum (ascending).
    cnt = cnt(I,:);
else
    cnt = opts.cnt;
end
%

% Label states
if opts.labelstates
    disp(['Labeling states for ' num2str(numel(cc)) ' runs...']);
    for i = 1:numel(cc)
        d{i} = pdist2(tanh(cc{i})',cnt,'correlation'); % label each full matrix using output centroids
        [dmin{i},st{i}] = min(d{i},[],2); % state is determined by minimum dist to centroids.
    end
end

% Determine state probabilities
% if opts.calc_prob
%     disp(['Getting state probabilities for ' num2str(numel(cc)) ' runs...']);
%     p = [];
%     for i = 1:numel(cc)
%         for j = 1:size(cc{i},2)
%             V1 = cc{i}(:,j);
%             V1 = V1/norm(V1);
%             for k = 1:opts.k
%                 V2 = cnt(k,:)';
%                 V2 = V2/norm(V2);
%                 p{i}(j,k) = (V1' * V2)^2;
%             end
%         end
%     end
% end
disp('Done')