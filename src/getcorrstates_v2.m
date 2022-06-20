
function [results, opts] = getcorrstates_v2(runnames,opts)
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

% rng reset for consistency
rng(555);
% Set default options
showprogbar = 1;
if ~isfield(opts,'labelstates'); opts.labelstates = 1; end
if ~isfield(opts,'nreps'); opts.nreps = 1; end
if isstr(runnames); runnames = {runnames};  showprogbar = 0; end

results = []; cc_ds = [];
warning('off','all')

% Compose correlation matrices
opts.win = getwin(opts.ww,opts.sig); % Create window (square if sig = 0, tapered otherwise)
if showprogbar; h = waitbar(0,['Getting sliding correlations for ' num2str(numel(runnames)) ' runs...']); h.Position(1:2) = [0 1080]; end
for i = 1:numel(runnames) % iterate through each cell of H.
    if ~strcmp(runnames{i}(end-5:end),'_H.mat')
        runname = [runnames{i} '_H.mat'];
    else
        runname = runnames{i};
    end
    load(fullfile(opts.Hdir,runname),'H','m') % Load data
    if opts.clustervar == 'j'
        Hall{i} = H.jrgeco; 
    else
        Hall{i} = H.chbo+H.chbr; 
    end
    rotf{i} = m.rotf;
    cc_ds{i} = slcorr(Hall{i},opts.skp,opts.win,rotf{i},1,0);
    if showprogbar; waitbar(i/numel(runnames),h); end
end
if showprogbar; close(h); end
cc_ds = cell2mat(cc_ds)';
%

% Perform clustering
if opts.dokmeans
    if ~isfield(opts,'cnt')
        [st,opts.cnt] = kmeans(cc_ds,opts.k,'Distance','sqeuclidean','Replicates',opts.nreps);
    else
        [st,opts.cnt] = kmeans(cc_ds,opts.k,'Distance','sqeuclidean','Start',opts.cnt);
    end

    for i = 1:opts.k
       opts.cnt(i,:) = mean(cc_ds(find(st==i),:)); 
    end
end
%

% Label states
if opts.labelstates
    if showprogbar; h = waitbar(0,['Labeling states for ' num2str(numel(runnames)) ' runs...']); h.Position(1:2) = [0 1080]; end
    for i = 1:numel(runnames)
        [cc,results{i}.idx] = slcorr(Hall{i},1,opts.win,rotf{i},1,0);
        if isempty(cc)
            continue
        end
        results{i}.d = [];
        results{i}.st = [];
        results{i}.Mr = [];
        for r = 1:size(cc,2)
            dtmp = lsqnonneg(opts.cnt',cc(:,r));
            [~,sttmp] = max(dtmp);
            rtmp = corr(opts.cnt(sttmp,:)',cc(:,r));
            results{i}.d(r,:) = dtmp;
            results{i}.st(r) = sttmp;
            results{i}.Mr(r) = 1-rtmp;
        end
        results{i}.rotf = rotf{i};
        if showprogbar; waitbar(i/numel(runnames),h); end
    end
    if showprogbar; close(h); end
end

