function [H,IDXout,nmerge] = merge_H(H,IDX,thr)
nmerge = 0;
nmerge_prev = NaN;
while 1
    cc = corrcoef(H'); cc = cc.^2;
    if max(cc(:)) < thr || nmerge_prev == nmerge
        break
    else
        nmerge_prev = nmerge;
        cc(cc < thr) = 0;
        cc(cc==1) = 0;
        done = 0;
        while ~done
            [max_r,max_c] = find(cc==max(cc(:)));
            maxccIDX = double(IDX==max_r(1)) + double(IDX==max_r(2));
            bwcc = bwconncomp(maxccIDX);
            if bwcc.NumObjects == 1
                r1 = max_r(1); r2 = max_r(2);
                s1 = numel(find(IDX==r1));
                s2 = numel(find(IDX==r2));
                H(r1,:) = (H(r1,:)*s1 + H(r2,:)*s2)/(s1+s2);
                H(r2,:) = 0;
                IDX(IDX==r2) = r1;
                nmerge = nmerge+1;
                done = 1;
            else
                cc(max_r,max_c) = 0;
            end
            if nansum(cc(:)) == 0
                done = 1;
            end
        end
    end
end
un = unique(IDX(:)); un(un==0) = [];
IDXout = IDX;
for i = 1:numel(un)
   IDXout(IDX==un(i)) = i;
end
H = H(un,:);
