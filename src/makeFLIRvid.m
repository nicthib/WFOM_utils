% [saved,opts] = makeFLIRvid(opts) converts a stream of jpegs to an avi,
% with an added crop functionality.
% 
% Required inputs (fields of the input struct 'opts'):
% mouse: String indicating the mouse id and day. (example: 'cm124_1')
% run: string indicating the run you want to process (example: 'runB')
% tag: appended name to filename (example: 'pupil')
% frames: 2 element vector with start and end frame. if empty, loads the
% whole run
% savepath: leave empty to save to dingus/1/RS_analysis/webcams
%
% Outputs:
% saved: full path and name of saved video. Use the command
% playvideo(saved) to view the video.
% opts: The modified opts struct 

function [saved,opts] = makeFLIRvid(opts)
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
if ~isfield(opts,'savepath')
    opts.savepath  = '/local_mount/space/dingus/1/RS_analysis/webcams';
end
if ~isfield(opts,'frames')
    framecommand = [];
    opts.frames = 0;
else
    framecommand = ['-frames:v ' num2str(abs(diff(opts.frames)))];
end
if ~isfield(opts,'tag')
    opts.tag = '';
end
if ~isfield(opts,'framerate')
    opts.framerate = 59.9; % For 2019 RS jrgeco cohort
end
if ~isfield(opts,'timestamp')
    opts.timestamp = 0;
end
if ~isfield(opts,'format')
    opts.format = '.avi';
end
if ~isfield(opts,'mouse')
    disp('Please input a mouse id into opts.mouse!')
    return
end
if ~isfield(opts,'run')
    disp('Please input a run name into opts.run!')
    return
end

% IMPORTANT STEP: This adds the static ffmpeg build in home/nic to
% the PATH environment variable so you can access it in the system command.
% See dingus/1/RS_analysis/ffmpeg-readme.txt file for instructions on 
% how to install your own static build.
buildname = ':/home/evan/ffmpeg-git-20200119-amd64-static';
if isempty(regexp(getenv('PATH'),buildname)) 
    setenv('PATH', [':' getenv('PATH') buildname]);
end

nc = '10'; % number of cores to use in ffmpeg
% im0dsf = dsf for picking crop. This is in the edge case that the two 
% webcams are capturing at different resolutions.
im0dsf = 1; im1dsf = 1; cams = [];
loadpath = fullfile(findmousefolder(opts.mouse),'webcam',opts.run);
% Set save location/name
saved = fullfile(opts.savepath,[opts.mouse '_' opts.run '_' opts.tag opts.format]);
% Picking crop area
try % Tries to load cam0. Need to improve this to avoid try statement
    imtmp = LoadFLIR(loadpath,1,0); im0 = imresize(imtmp,[540 720]);
    im0dsf = size(imtmp,1)/size(im0,1);
    cams = [cams 0];
    % ^This is to preserve original frame size for eventual image cropping^
catch
    im0 = [];
end
try % Tries to load cam1. Need to improve this to avoid try statement
    imtmp = LoadFLIR(loadpath,1,1); im1 = imresize(imtmp,[540 720]);
    im1dsf = size(imtmp,1)/size(im1,1);
    cams = [cams 1];
catch
    im1 = [];
end
if ~isfield(opts,'x')
    figure('Position',[0 0 1920 1080]/2)
    imagesc(cat(2,im0,im1)); axis image; axis off
    title('Click twice to crop. Don''t overlap cams!')
    [x,y] = ginput(2);
    opts.x = x; opts.y = y;
else
    x = opts.x; y = opts.y;
end
% Edge case of two webcams, and you picked the right one
if x(1) > 720 && numel(cams) == 2
    cam = 1; x = x-720;
    % Two webcams, but you picked the left one
elseif x(1) <= 720 && numel(cams) == 2
    cam = 0;
    % Otherwise, the cam must be the only one loaded into the cams variable.
else
    cam = cams;
end
% Rescale coordinates to full resolution of original image
if cam == 0
    x = x*im0dsf; y = y*im0dsf;
elseif cam == 1
    x = x*im1dsf; y = y*im1dsf;
end
x = round(x); y = round(y);
vx = num2str(min(x)); w = num2str(abs(diff(x))); % (vx,vy) is top left corner of crop
vy = num2str(min(y)); h = num2str(abs(diff(y))); % (h,w) is height and width
close all

if opts.timestamp
    ts = [', drawtext=\timecode=''' maketimestamp(opts.frames(1),60) ''':rate=60:fontcolor=white'];
else
    ts = '';
end

if isfield(opts,'height')
    heightcommand=[' , scale=-1:' num2str(opts.height) ',crop=trunc(iw/2)*2:trunc(ih/2)*2'];
else
    heightcommand='';
end

% ffmpeg flags used:
% -y: force overwrite
% -framerate: obvious
% -start_number: image number to start at
% -threads: cpu threads to use when writing
% -i: input file (%d is the image index number)
% -frames:v: number of frames to write to file (framecommand variable)
% -vcodec: video codec (mjpeg makes it vieable on the server)
% -vf: video filter (I used crop=, drawtext=, and \timecode=
% -q:v: video quality (1 is best?)
ffmpeg_command = ['ffmpeg -loglevel panic -y -framerate ' mat2str(opts.framerate) ' -start_number ' num2str(opts.frames(1)) ' -threads ' nc ' -i ' loadpath '/' opts.run '_stim1_%d_cam'...
    mat2str(cam) '.jpg ' framecommand ' -vcodec mjpeg -vf "crop=' w ':' h ':' vx ':' vy ...
    ' ' ts heightcommand ' " -b 20M -loglevel repeat ' saved];
system(ffmpeg_command)
