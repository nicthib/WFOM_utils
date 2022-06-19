function [L,R] = midline(pts)
midptx = pts(:,1); midpty = pts(:,2);
if numel(midptx) ~= 2
    errordlg('Please pick only 2 points!')
    return
end

% re-evaluate line at top and bottom of frame
ri = diff(midpty);
ru = diff(midptx);
sl = ri/ru;
b = midpty(1)-sl*midptx(1);
midpty = [0 256];
midptx = ([0 256] - b)/sl;
L = poly2mask([midptx 0 0],[midpty 256 0],256,256);
R = ~L;