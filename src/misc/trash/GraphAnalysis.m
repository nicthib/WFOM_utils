H_tr = H_nf{3};
H_tr_z = [];
zcl = 11980;
for i = 0:17
    H_tr_z(:,i*zcl+1:(i+1)*zcl) = zscore(H_tr(:,i*zcl+1:(i+1)*zcl)')';
end
for i = 1:size(H_tr_z)
    if any(isnan(H_tr_z(i,:)))
        H_tr_z(i,:) = NaN;
    end
end
H_tr_z(find(sum(H_tr_z,2) == 0),:) = [];
H_tr_z(find(isnan(H_tr_z(:,1))),:) = [];
H_tr_z=sepblockfun(H_tr_z,[1,5],@mean);
%% Threshold

H_tr_z(H_tr_z > 1) = 1;
H_tr_z(H_tr_z < -1) = -1;
H_tr_z(H_tr_z < 1 & H_tr_z > -1) = 0;
A = H_tr_z';

%%
[C,~,IC] = unique(A,'rows');
% IC is the 
s = IC;
t = [IC(2:end)' 1];
G = digraph(s,t,ones(numel(s),1));
plot(G,'xdata',sum(a,2),'ydata')


G = simplify(G,'sum');
G = rmedge(G,find(G.Edges.Weight == 1));
G = rmnode(G,find(outdegree(G)+indegree(G) == 0));
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)