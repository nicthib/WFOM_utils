function [ord st_out] = testoverlap(st1,st2,k)
p = perms(1:k);
for i = 1:size(p,1)
    st_tmp = zeros(size(st2));
    for ks = 1:k
        st_tmp(st2==ks) = p(i,ks);
    end
    ov(i) = numel(find(st1-st_tmp==0))/numel(st1);
end
ord = p(min(find(ov==max(ov))),:);
disp(['Overlap score: ' mat2str(max(ov))])
for ks = 1:k
    st_out(st2==ks) = ord(ks);
end