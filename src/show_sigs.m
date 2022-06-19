function pairs = show_sigs(p)
pairs = [];
p = p+~triu(p,1);
reps = sum(sum(triu(ones(size(p)),1)));
for i = 1:size(p,1)
    for j = 1:size(p,2)
        if p(i,j) < .05/reps
            pairs(end+1,:) = [i j];
        end
    end
end
