function out = LoadFLIR(path,i,cam)
% Loads a single FLIR image from path, using i as frame # and cam as cam id (0 or 1).
tmp = dir(fullfile(path,['*_' mat2str(i) '_cam' mat2str(cam) '.jpg']));
if ~isempty(tmp)
    out = uint8(imread(fullfile(path,tmp.name)));
else
    disp(['FILE MISSING for frame # ' mat2str(i)])
    out = [];
end

