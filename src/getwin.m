function win = getwin(w,sig)
% win = getwin(w,sig) returns a window with width w and convolved with a gaussian kernel.
x = -1000:1000;
if sig > 0
    y = gaussmf(x,[sig 0]);
    y(y<.0001) = [];
    win = conv(ones(1,w),y);
    win = win/max(win);
else
    win = ones(1,w);
end