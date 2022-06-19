runname = 'runD';
a = findmousefolder('cm150_1');
wfiles = dir(fullfile(a,'CCD',runname,'*.jpg'));
parfor i = 1:numel(wfiles)
    movefile(fullfile(wfiles(i).folder,wfiles(i).name),fullfile(a,'webcam',runname))
end
disp('DONE')