% ts = maketimestamp(st_frame,fr) creates the timestamp for use in ffmpeg
% commands. standard format of timestamp is 00\:00\:00\:00. This command
% allows for non-zero starting timestamps when converting videos using
% ffmpeg.
%
% Inputs:
% st_frame: integer of starting frame #
% fr: frame rate (typically 60)
function ts = maketimestamp(st_frame,fr)
st_full = st_frame/fr;
st_m = num2str(floor(st_full/fr));
if numel(st_m) == 1
    st_m = ['0' st_m];
end
st_s = num2str(floor(mod(st_full,fr)));
if numel(st_s) == 1
    st_s =  ['0' st_s];
end
st_ss = num2str(round((st_full-floor(st_full))*100));
if numel(st_ss) == 1
    st_ss =  ['0' st_ss];
end
ts = ['00\:' st_m '\:' st_s '\:' st_ss];
