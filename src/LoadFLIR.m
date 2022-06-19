function out = LoadFLIR(path,i,cam)
tmp = dir(fullfile(path,['*_' mat2str(i) '_cam' mat2str(cam) '.jpg']));
if ~isempty(tmp)
    out = uint8(imread(fullfile(path,tmp.name)));
else
    disp(['FILE MISSING for frame # ' mat2str(i)])
    out = [];
end

