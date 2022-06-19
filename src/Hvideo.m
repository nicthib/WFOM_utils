%% Video from H and kmeans
n = 1;
IDX = cm125.IDXs_b{n};
H = H_n{n};
for i = 1:1000
    tmp = zeros(size(IDX));
    for j = 1:max(IDX(:))
        tmp(IDX == j) = H(j,i);
    end
    imagesc(tmp)
    axis image; axis off; caxis([-.1 .1])
    pause(.05)
end