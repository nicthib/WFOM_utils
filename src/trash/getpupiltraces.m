function [xc,yc,R] = getpupiltraces(data,thresh)
f = fieldnames(data);
j = 1;
for i = 1:numel(f)-1
    tmp = data.(f{i});
    if isa(tmp,'double') && size(tmp,2) == 3
        pc(:,j,1:2) = tmp(:,1:2);
        pc(find(tmp(:,3) < thresh),j,:) = NaN;
        j = j+1;
    end
end

for i = 1:8
    for j = 1:2
        pcs(:,i,j) = sepblockfun(squeeze(pc(1:end-mod(size(pc,1),3),i,j)),3,'nanmean');
    end
end

for i = 1:size(pcs,1)
    [xc(i),yc(i),R(i)] = circfit(squeeze(pcs(i,:,1)),squeeze(pcs(i,:,2)));
end