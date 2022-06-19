function delrows = cleanupH(H,delrows)
H(isnan(H)) = 0;
H(isinf(H)) = 0;
for i = 1:size(H,1)
    if sum(H(i,:)) == 0
        delrows = [delrows i];
    end
end

