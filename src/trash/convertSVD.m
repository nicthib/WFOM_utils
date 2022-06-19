function out = convertSVD(C,S,comps,nanidx,sz)
tmp = zeros(sz^2,size(S,2));
tmp(nanidx,:) = C(:,comps)*S(comps,:);
tmpmin = min(tmp(:));
if tmpmin < 0
    tmp = tmp-tmpmin;
end
tmpstd = find(std(tmp(:,30:100),[],2)<1e-10);
tmp(tmpstd,:) = NaN;
out = reshape(tmp,[sz sz size(S,2)]);