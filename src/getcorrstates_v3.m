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

function [results, opts] = getcorrstates_v3(runnames,opts)
rng(555);
% Set default options
showprogbar = 1;
if ~isfield(opts,'labelstates'); opts.labelstates = 1; end
if ~isfield(opts,'nreps'); opts.nreps = 1; end
if isstr(runnames); runnames = {runnames};  showprogbar = 0; end
results = cell(numel(runnames),1);

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
    Hj{i} = H.jrgeco;
    Hh{i} = H.chbo+H.chbr; 
    do_detrend = 1;
    rotf{i} = m.rotf(:);
    [ccj{i},results{i}.idx] = slcorr(Hj{i},opts.skp,opts.win,rotf{i},do_detrend);
    cch{i} = slcorr(Hh{i},opts.skp,opts.win,rotf{i},do_detrend);
    if showprogbar; waitbar(i/numel(runnames),h); end
end
if showprogbar; close(h); end
ccj = cell2mat(ccj)';
cch = cell2mat(cch)';
%

% Perform clustering
if opts.dokmeans
    if ~isfield(opts,'cntj')
        [opts.stj,opts.cntj] = kmeans(ccj,opts.k,'Distance','sqeuclidean','Replicates',opts.nreps);
    else
        [opts.stj,opts.cntj] = kmeans(ccj,opts.k,'Distance','sqeuclidean','Start',opts.cntj);
    end

    % hemo states
    if opts.doh
        for i = 1:opts.k
            opts.cnth(i,:) = mean(cch(find(opts.stj==i),:)); 
        end
    end
end
%

% Label states
if opts.labelstates
    if showprogbar; h = waitbar(0,['Labeling states for ' num2str(numel(runnames)) ' runs...']); h.Position(1:2) = [0 1080]; end
    for i = 1:numel(runnames)
        [ccj,results{i}.idx] = slcorr(Hj{i},1,opts.win,rotf{i},do_detrend);
        if isempty(ccj)
            continue
        end
        n = size(ccj,2);
        results{i}.dj = zeros(n,opts.k);
        results{i}.stj = zeros(n,1);
        for r = 1:n
            if opts.nnls
                dtmp = lsqnonneg(opts.cntj',ccj(:,r));
                [~,sttmp] = max(dtmp);
            else
                dtmp = pdist2(opts.cntj,ccj(:,r)');
                [~,sttmp] = min(dtmp);
            end
            results{i}.dj(r,:) = dtmp;
            results{i}.stj(r) = sttmp;
        end
        
        if opts.doh
            [cch,results{i}.idx] = slcorr(Hh{i},1,opts.win,rotf{i},do_detrend);
            n=size(cch,2);
            results{i}.dh = zeros(n,opts.k);
            results{i}.sth = zeros(n,1);
            for r = 1:n
                if opts.nnls
                    dtmp = lsqnonneg(opts.cnth',cch(:,r));
                    [~,sttmp] = max(dtmp);
                else
                   dtmp = pdist2(opts.cnth,cch(:,r)');
                    [~,sttmp] = min(dtmp);
                end
                results{i}.dh(r,:) = dtmp;
                results{i}.sth(r) = sttmp;
            end
        end
        results{i}.rotf = rotf{i};
        if showprogbar; waitbar(i/numel(runnames),h); end
    end
    if showprogbar; close(h); end
end

