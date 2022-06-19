function show_state_centroid(varargin)
cc = varargin{1};
r = sqrt(numel(cc(:)));
cc = reshape(cc,[r r]);

% defaults
cIDX_comps = [];
ticklabels = '';
cmap_labels = lines(numel(cIDX_comps));

if numel(varargin) > 1
    cIDX_comps = varargin{2};
end

if ~isempty(cc)
    
end
caxis([0 1])
% for i = 1:size(cmap_labels,1)
%     tmp = mat2str(cmap_labels(i,:));
%     tmp = tmp(2:end-1);
%     ticklabels{i} = ['\color[rgb]{' tmp '} ' ticklabels{i}];
% end
% xticks(laboffset);
% yticks(laboffset);
% xticklabels(ticklabels);
% yticklabels(ticklabels);
% xtickangle(45)

axis image; colormap(jet(128))

xlim([.5 r+.5]); ylim([.5 r+.5])
