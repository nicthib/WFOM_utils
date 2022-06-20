% assign_IDX() assigns values in val (vector with N values) for IDX redions 1 to N.
function IDXout = assign_IDX(IDX,val)
IDXout = zeros(size(IDX));
for i = 1:numel(val)
    IDXout(IDX==i) = val(i);
end