function highlight_network(IDX,network,color)
network_BW = ismember(IDX,network);
contour(network_BW,1,'Color',color,'LineWidth',2);