function t = d_to_t(d1,d2)
t = zeros(numel(d1));
for i = 1:numel(d1)
    for j = 1:numel(d2)
        t(i,j) = d2(i)-d1(j);
    end
end
