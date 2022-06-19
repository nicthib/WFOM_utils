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

function [results, opts] = getcorrstates_behav(runnames,opts)

% Set default options
if ~isfield(opts,'variable'); opts.variable = 'jrgeco'; end
if ~isfield(opts,'nreps'); opts.nreps = 1; end
if ~isfield(opts,'Hdir'); opts.Hdir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; end
if ~isfield(opts,'highpassH'); opts.highpassH = 1; end
if isstr(runnames); runnames = {runnames}; end
%

% Inits
results = [];
cnt_j_full = cell(5,1); cnt_h_full = cell(5,1);
warning('off','all') % K-means tends to send warnings while in parallel.
opts.win = getwin(opts.ww,opts.sig); % Create window (square if sig = 0, tapered otherwise)
opts.k = 5; % state k must be 5 (for now)
%

% Load data, get states and compose exemplar FC maps
h = waitbar(0,['Getting sliding correlations for ' num2str(numel(runnames)) ' runs...']); h.Position(1:2) = [0 1080];
for i = 1:numel(runnames) % iterate through each cell of H.
    % Load, assign
    load(fullfile(opts.Hdir,runnames{i}),'H','m') % Load data
    Hj = H.jrgeco;
    Hh = H.chbo+H.chbr;
    %
    
    % Filtercnt_j_full 
    bp = 20/(opts.wr*2);
    if opts.highpassH
        Hj = highpass(Hj',bp,20)';
        Hh = highpass(Hh',bp,20)';
    end
    Hh = lowpass(Hh',2,20)';
    %
    
    % Get states from rotf
    [st,st_idx,rotf_st] = behav_states(m.rotf,opts.ww,1);
    %
    
    % Get FCs and append to full state matrix
    for s = 1:numel(st)
        idx = st_idx(s)-opts.wr:st_idx(s)+opts.wr;
        cnt_j_full{st(s)}(end+1,:) = colvec(corrcoef(Hj(:,idx)'));
        cnt_h_full{st(s)}(end+1,:) = colvec(corrcoef(Hh(:,idx)'));
    end
    waitbar(i/numel(runnames),h)
    %
end
close(h)
%

% Average all FC instances for each state
for k = 1:opts.k
    opts.cnt(k,:) = mean(cnt_j_full{k},1);
    opts.cnth(k,:) = mean(cnt_h_full{k},1);
end

%
% 
% % Label states
% if opts.labelstates
%     h = waitbar(0,['Labeling states for ' num2str(numel(runnames)) ' runs...']); h.Position(1:2) = [0 1080];
%     for i = 1:numel(runnames)
%         if opts.donnmf
%             for r = 1:size(cc_j{i},1)
%                 results{i}.dj(r,:) = lsqnonneg(opts.cnt',cc_j{i}(r,:)');
%                 results{i}.dh(r,:) = lsqnonneg(opts.cnth',cc_h{i}(r,:)');
%             end
%             [~,results{i}.st] = max(results{i}.dj,[],2);
%             
%             waitbar(i/numel(runnames),h)
%         else
%             results{i}.dj = pdist2(cc_j{i},opts.cnt,'correlation'); % label each full matrix using output centroids
%             results{i}.dh = pdist2(cc_h{i},opts.cnth,'correlation'); % label each full matrix using output centroids
%             [~,results{i}.st] = min(results{i}.dj,[],2); % state is determined by minimum dist to centroids.
%             waitbar(i/numel(runnames),h)
%         end
%     end
%     close(h)
% end

