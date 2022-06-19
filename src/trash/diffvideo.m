function out = diffvideo(file)
vidObj = VideoReader(file);
im = [];
while hasFrame(vidObj)
    vF = readFrame(vidObj);
    im(:,:,end+1) = vF(:,:,1);
end
imd = diff(im,1,3); imd = reshape(imd,[size(imd,1)*size(imd,2) size(imd,3)]);
out = std(single(imd),[],1);