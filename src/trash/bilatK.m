function nb = bilatK(IDX)
IDX = uint16(IDX);
for i = 1:max(IDX(:))
    tmp = IDX == i;
    bw = bwconncomp(bwareaopen(tmp,50));
    if numel(bw.PixelIdxList) > 1
        nb(i) = 1;
    else
        nb(i) = 0;
    end
end