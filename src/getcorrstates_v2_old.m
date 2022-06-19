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

function [results, opts] = getcorrstates_v2(runnames,opts)

% Set default options
showprogbar = 1;
if ~isfield(opts,'labelstates'); opts.labelstates = 1; end
if ~isfield(opts,'nreps'); opts.nreps = 1; end
if ~isfield(opts,'donnmf'); opts.donnmf = 0; end
if isstr(runnames); runnames = {runnames};  showprogbar = 0; end

results = []; cc_dsj = [];
warning('off','all') % K-means tends to send warnings while in parallel.

% Compose correlation matrices
opts.win = getwin(opts.ww,opts.sig); % Create window (square if sig = 0, tapered otherwise)
if showprogbar; h = waitbar(0,['Getting sliding correlations for ' num2str(numel(runnames)) ' runs...']); h.Position(1:2) = [0 1080]; end
for i = 1:numel(runnames) % iterate through each cell of H.
    if ~strcmp(runnames{i}(end-5:end),'_H.mat')
        runname = [runnames{i} '_H.mat'];
    else
        runname = runnames{i};
    end
    load(runname,'H','m') % Load data
    Hh{i} = H.chbt; Hj{i} = H.jrgeco; rotf{i} = m.rotf;
    cc_dsj{i} = slcorr(Hj{i},opts.skp,opts.win,rotf{i});
    cc_dsh{i} = slcorr(Hh{i},opts.skp,opts.win,rotf{i});
    if showprogbar; waitbar(i/numel(runnames),h); end
end
if showprogbar; close(h); end
cc_dsj = cell2mat(cc_dsj)';
cc_dsh = cell2mat(cc_dsh)';
%

% Perform clustering
options = statset('UseParallel',0); % To increase speed
distfun = 'sqeuclidean';
if strcmp(distfun,'correlation')
    distfun2 = 'correlation';
else
    distfun2 = 'squaredeuclidean';
end
if opts.dokmeans
    if opts.clusterj
        if ~isfield(opts,'cntj')
            [stj,opts.cntj] = kmeans(cc_dsj,opts.k,'Distance',distfun,'options',options,'Replicates',opts.nreps);
        else
            [stj,opts.cntj] = kmeans(cc_dsj,opts.k,'Distance',distfun,'options',options,'Start',opts.cntj);
        end

        for i = 1:opts.k
           opts.cnth(i,:) = mean(cc_dsh(find(stj==i),:)); 
        end
    else
        if ~isfield(opts,'cnth')
            [sth,opts.cnth] = kmeans(cc_dsh,opts.k,'Distance',distfun,'options',options,'Replicates',opts.nreps);
        else
            [sth,opts.cnth] = kmeans(cc_dsh,opts.k,'Distance',distfun,'options',options,'Start',opts.cnth);
        end

        for i = 1:opts.k
           opts.cntj(i,:) = mean(cc_dsj(find(sth==i),:)); 
        end
    end
end
%

% Label states
if opts.labelstates
    if showprogbar; h = waitbar(0,['Labeling states for ' num2str(numel(runnames)) ' runs...']); h.Position(1:2) = [0 1080]; end
    for i = 1:numel(runnames)
        [cc_j,results{i}.idx] = slcorr(Hj{i},1,opts.win,rotf{i},1);
        cc_h = slcorr(Hh{i},1,opts.win,rotf{i},1);
        if isempty(cc_j)
            continue
        end
        realidx = find(~isnan(cc_j(:,1)));
        if opts.donnmf
            for r = 1:size(cc_j,2)
                results{i}.dj(r,:) = lsqnonneg(opts.cntj(:,realidx)',cc_j(realidx,r));
                results{i}.dh(r,:) = lsqnonneg(opts.cnth(:,realidx)',cc_h(realidx,r));
            end
            [~,results{i}.stj] = max(results{i}.dj,[],2);
            [~,results{i}.sth] = max(results{i}.dh,[],2);
        else
            results{i}.dj = pdist2(cc_j(realidx,:)',opts.cntj(:,realidx),distfun2);
            results{i}.dh = pdist2(cc_h(realidx,:)',opts.cnth(:,realidx),distfun2);
            [~,results{i}.stj] = min(results{i}.dj,[],2); 
            [~,results{i}.sth] = min(results{i}.dh,[],2); 
        end
        if showprogbar; waitbar(i/numel(runnames),h); end
    end
    if showprogbar; close(h); end
end

