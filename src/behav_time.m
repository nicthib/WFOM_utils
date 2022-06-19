function [ttn,ttp] = behav_time(TC)
TC = TC(:);
% TC is logical valuie of behavior
TCf = flipud(TC);
p = 0; n = 0;
% previous
for i = 1:numel(TC)
    if TC(i)
        p = 0;
    else
        p = p+1;
    end
    ttp(i) = p;
    
    if TCf(i)
        n = 0;
    else
        n = n+1;
    end
    ttn(i) = n;
end
ttp = ttp(:);
ttn = flipud(ttn(:));