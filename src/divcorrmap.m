function c = divcorrmap(c,divs)
if ~isempty(divs)
    r = size(c,1);
    n = 3;
    c = [c(:,1:divs(1)) nan(r,n) c(:,divs(1)+1:divs(2)) nan(r,n) c(:,divs(2)+1:end)];
    c = [c(1:divs(1),:); nan(n,r+n*2); c(divs(1)+1:divs(2),:); nan(n,r+n*2); c(divs(2)+1:end,:)];
end