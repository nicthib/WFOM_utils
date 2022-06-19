function label_corrmap(networks,ticklabels)
r = sum(networks);
regions_all = [cumsum(networks) cumsum(networks)+r];
laboffset = regions_all - [networks networks]/2 + .5;
cmap_labels = lines(numel(networks));
cmap_labels = [cmap_labels;cmap_labels];
for i = 1:size(cmap_labels,1)
    tmp = mat2str(cmap_labels(i,:));
    tmp = tmp(2:end-1);
    ticklabels{i} = ['\color[rgb]{' tmp '} ' ticklabels{i}];
end
xticks(laboffset);
yticks(laboffset);
xticklabels(ticklabels);
yticklabels(ticklabels);
xtickangle(45);