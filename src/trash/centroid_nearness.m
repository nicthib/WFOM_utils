function score = centroid_nearness(IDX1,IDX2)

for i = 1:max(IDX1(:))
    tmp1 = IDX1==i;
    tmp2 = IDX2==i;
    s1 = regionprops(tmp1,'centroid');
    s2 = regionprops(tmp2,'centroid');
    c1 = s1.Centroid;
    c2 = s2.Centroid;
    sc(i) = sqrt(sum((c1-c2).^2));
end
score = mean(sc);