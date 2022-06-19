% bl = baselinefromrot(rot,span,thr) obtains a baseline from rotary data
% 'rot', defined as at least 'span' values below 'thr'. In essence, this
% function looks for some frame of time where the rotary value is not
% changing for at least 'span' frames. the returned variable 'bl'
% represents the frames where a baseline was found.
function bl = baselinefromrot(rot,span,thr)
rot = rot(:); % vectorize
for i = 1:numel(rot)-span-1
    idx = i:i+span-1;
    if max(rot(idx)) < thr && ~any(isnan(rot(idx)))
        bl = idx;
        return
    end
end
bl = [];
disp('No baseline found!')
     