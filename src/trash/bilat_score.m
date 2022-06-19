function b_score = bilat_score(cc)
half_sz = size(cc,1)/2;
idx = 1:half_sz;
for i = 1:numel(idx)
    b_score(i) = cc(idx(i),idx(i)+half_sz);
end