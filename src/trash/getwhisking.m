function s = getwhisking(mouse,run,check,s)
try parpool(6); catch; end
path = findmousefolder(mouse);
path = fullfile(path,'webcam',run);
nframes = (numel(dir(path))-2);
cam = 0;
if check == 1
    im = LoadFLIR(path,1,cam);
    imagesc(im);
    [x,y] = ginput(2);
    x = uint16(x); y = uint16(y);
else
    x = s.x; y = s.y;
end
close all
im = zeros(diff(y)+1,diff(x)+1,nframes,'int16');
tic
ppm = ParforProgMon('Getting whiskers...', nframes);
parfor i = 1:nframes
    tmp = LoadFLIR(path,i,cam);
    if ~isempty(tmp)
        im(:,:,i) = tmp(y(1):y(2),x(1):x(2),:);
    else
        im(:,:,i) = 0;
    end
    ppm.increment();
end
fprintf('Done, fr = %.1f\n',nframes/toc)

imd = diff(im,1,3); imd = reshape(imd,[201^2 size(imd,3)]);
s.whisk = std(single(imd),[],1);
s.im = im(:,:,1:10);
s.x = x; s.y = y;

