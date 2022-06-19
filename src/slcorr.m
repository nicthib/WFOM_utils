% cc = slcorr(H,skp,win) calculates a sliding correlation matrix.
% 
% Inputs: 
% H: an n x t matrix, where n is the number of spatial components and t isd
% the number of frames.
% 
% skp: the skip value (this is equivalent to opts.skipfactor in
% getcorrstates())
% 
% win: the sliding window to be used (obtained with the function getwin())
%
% Outputs:
% 
% cc: the correlation matrix output. Note that it will decrease in size as
% a function of increasing the skip factor.

function [cc,cluster_frames] = slcorr(varargin)
H = varargin{1};
skp = varargin{2};
win = varargin{3};
rotf = varargin{4};
do_detrend = varargin{5};
delay = varargin{6};

% pad array
ww = numel(win);
wr = (ww-1)/2;
r = size(H,1);
n = size(H,2);
cc = [];
cluster_bool = zeros(n,1);
cluster_frames = [];
W = 1000; % 30 seconds + radius

% normalize H
%H = zscore(H')';

% First, get cluster frames
if isempty(rotf)
    return
elseif strcmp(rotf,'all')
    cluster_bool = ones(1,size(H,2));
else
    epochs = extract_session(rotf);
    if numel(epochs) < 1
        return
    end
    for i = 1:size(epochs,1)
        if (epochs(i,2) - epochs(i,1)) >= W
            cluster_bool((epochs(i,1)+ww):skp:(epochs(i,2)-ww-1)) = 1;
        end
    end
end
cluster_frames = find(cluster_bool==1);
cluster_frames(cluster_frames<(wr+1)) = [];
cluster_frames(cluster_frames>(size(H,2)-wr)) = [];

if mod(ww,2) == 1
    cc = zeros(r,r,numel(cluster_frames));
    fwin = repmat(win,[r 1]); % repmat so that win can be applied to every value in H
    for i = 1:numel(cluster_frames)
        idx = [-wr:wr] + cluster_frames(i) + delay;
        tmp = H(:,idx).*fwin; % Multiply by tapered window
        if do_detrend
            tmp = detrend(tmp'); % Detrend
        else
            tmp = tmp';
        end
        cc(:,:,i) = corr(tmp);
    end
    cc = reshape(cc,[r^2 size(cc,3)]);
else
    disp('Window must be an odd number width')
    return
end
