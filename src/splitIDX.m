function [L,Rspan,C,IDXout] = splitIDX(IDX)
%  Light GUI for simple splitting of ROI map. Obsolete in current workflow.
close all; figure;
nanidx = IDX == 0;
IDXtmp = zeros(size(IDX));
IDX = cleanupIDX(IDX,100); axis image
L = []; R = []; C = []; B = []; n = 1;
for i = 1:max(IDX(:))
    CC = bwconncomp(IDX==i);
    for j = 1:numel(CC.PixelIdxList)
        IDXtmp(CC.PixelIdxList{j}) = n;
        cla; imagesc(logical(IDXtmp == n) + logical(imgradient(IDX)));
        title([mat2str(i) '/' mat2str(max(IDX(:))) ' - Up = Center, Down = Bilateral'])
        axis image; axis off
        k = waitforbuttonpress;
        v = double(get(gcf,'CurrentCharacter'));
        switch v
            case 28 % L arrow
                L(end+1) = n;
            case 29 % R arrow
                R(end+1) = n;
            case 30 % U arrow
                C(end+1) = n;
            case 31 % D arrow
                B(end+1) = n;
            otherwise
                IDXtmp(CC.PixelIdxList{j}) = 0;
                n = n-1;
        end
        n = n+1;
    end
end
close all
IDXL = zeros(size(IDX));
IDXR = zeros(size(IDX));
IDXC = zeros(size(IDX));

n = 1;
for i = 1:numel(L)
    IDXL(IDXtmp == L(i)) = n;
    n = n+1;
end
for i = 1:numel(C)
    IDXC(IDXtmp == C(i)) = n;
    n = n+1;
end
for i = 1:numel(R)
    IDXR(IDXtmp == R(i)) = n;
    n = n+1;
end
try IDXL = sortkmeans(IDXL); catch; end
try IDXR = sortkmeans(IDXR); catch; end
try IDXC = sortkmeans(IDXC); catch; end
IDXout = IDXC+IDXL+IDXR;
IDXout(nanidx) = NaN;

%This code may have been artificially reducing corr values
%[zr,zc] = find(IDXout == 0);
%[nzr, nzc] = find(IDXout > 1);

%for i = 1:numel(zr)
%    [~,idx] = min(sum(([nzr,nzc]-repmat([zr(i),zc(i)],[numel(nzr),1])).^2,2));
%    IDXout(zr(i),zc(i)) = IDXout(nzr(idx),nzc(idx));
%end
