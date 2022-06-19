function IDXout = corrIDXmap(IDX,c)
IDXout = zeros(size(IDX));
for i = 1:max(IDX(:))
    IDXout(IDX==i) = mean(c(i,:));
end