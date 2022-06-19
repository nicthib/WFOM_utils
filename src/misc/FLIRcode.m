% Snippet of code to extract some behavioral epochs, and then look at
% them as a video
clear opts
opts.mouse = 'cm128_3';
opts.run = 'runB';
opts.savepath = '/local_mount/space/dingus/1/RS_analysis/webcams/epochs';
opts.timestamp = 1; % Add a timestamp to video
nframes = 60; % # of frames back and forward from epoch frame to make video.

% Some frames of interesting events. NOTE: If these frame numbers are from 
% the wfom data, they need to be multiplied by 3 to match the framerate 
% of the behavioral cam. Check the timesamp to make sure you're getting the
% right epoch time.
idx = [60 600 3000 12000]; 

for i = 1:numel(idx)
    opts.frames = [-nframes nframes]+idx(i); % Set frame range for video
    opts.tag = ['fs' num2str(i)]; % fs for "false starts" with the idx number for reference
    [~,opts] = makeFLIRvid(opts); % Return opts so that you don't have to crop every time.
end