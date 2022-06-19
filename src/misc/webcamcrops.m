%% Set environment variable
if isempty(regexp(getenv('PATH'),'/home/nic/ffmpeg/ffmpeg-git-20200119-amd64-static'))
    setenv('PATH', [getenv('PATH') ':/home/nic/ffmpeg/ffmpeg-git-20200119-amd64-static']);
end
%% PAWS
cd('S:/webcam/')
files = dir;
files = {files.name};

for i = 3:numel(files)
    
    if flip(i-2) == 0
        system(['ffmpeg -i ' files{i} ' -b:v 4M -filter:v "crop=720:270:1:270" ' files{i}(1:end-4) '_c.avi'])
    else
        system(['ffmpeg -i ' files{i} ' -b:v 4M -filter:v "crop=720:270:720:270" ' files{i}(1:end-4) '_c.avi'])
    end
end

%% PUPIL
savedir = 'S:/webcam/pupil';
cam = '1';
mouse = 'cm155_1';
runall = 'ABCDEFGHIJ';

for i = 1:numel(runall)
    run = ['run' runall(i)];
    loaddir = ['S:/' mouse '/webcam/' run];
    cd(loaddir)
    vidname = fullfile(savedir,[mouse '_' run '_p.avi']);
    if i == 2
        imagesc(imread(fullfile(loaddir,[run '_stim1_1_cam' cam '.jpg'])))
        axis image
        colormap gray
        [x,y] = ginput(1); x = x-75; y = y-75;
        close all
    end
    system(['ffmpeg -framerate 60 -start_number 0 -i ' run '_stim1_%d_cam' cam '.jpg -vcodec mpeg4 -q:v 1  -filter:v "crop=150:150:' mat2str(x) ':' mat2str(y) '" ' vidname])
end