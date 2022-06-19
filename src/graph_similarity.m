function s = graph_similarity(IDX1,IDX2)
denom = 0;
num = 0;
for i = 1:numel(IDX1)
    denom = denom + max(numel(IDX1{i}),numel(IDX2{i}));
    num = num + sum(ismember(IDX1{i},IDX2{i}));
end
s = num/denom;