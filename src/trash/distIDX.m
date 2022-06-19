function dIDX = distIDX(IDX,aiflag)
dIDX = zeros(256^2,1);
h = waitbar(0,'Getting distances...');
for i = 1:max(IDX(:))
    [X,Y] = meshgrid(1:256,1:256);
    b = [Y(:) X(:)];
    b(IDX == i,:) = [];
    a = reshape(IDX == i,[256 256]);
    [x,y] = ndgrid(1:size(a, 1), 1:size(a, 2));
    c = [x(logical(a)), y(logical(a))];
    if numel(c) > 2
        c = mean(c);
    end
    clear PQ
    [PQ(:,1), PQ(:,2)] = ind2sub([256 256], find(IDX == i));
    [~,bi] = dsearchn(b,PQ);
    if aiflag
        [~,ai] = dsearchn(c,PQ);
        dIDX(IDX==i) = (bi-ai)./max([ai bi],[],2);
    else
        dIDX(IDX==i) = bi;
    end
    waitbar(i/max(IDX(:)),h)
end
close(h)