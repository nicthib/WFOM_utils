function [l,y,g] = getstatelengths(st)
k = max(st);
l = cell(k,1);
y = []; g = [];
while numel(st) > 1
    if any(diff(st))
        tmp = min(find(diff(st) ~= 0)) + 1;
    else
        tmp = numel(st);
    end
    idx = st(1);
    l{idx}(end+1) = tmp - 1;
    st(1:tmp-1) = [];
    y(end+1) = tmp-1;
    g(end+1) = idx;
end

