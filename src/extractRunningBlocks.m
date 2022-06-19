function obj = extractRunningBlocks(obj)
% extractRunningBlocks separate the running bouts of each experimental session
% inputs: mousename, dataDIR (H_final directory), boutsDIR (a directory to
% save the output)
% output: runningBlocks (a structure array with fields like running onset, running offset,
% running duration, time to the previous running bout and time to the next running bout)
% input examples:
% mousename = "cm125";
% dataDIR = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final/';
% boutsDIR = '/home/bahar/Documents/Manuscript/results/runningBlocks/';

idx = 0;
% load one day, one run at a time
rz = getrunningpulses2(obj.B.rotf',0.01,20,100,"off");

% detect the running blocks (pulses)
rz_diff = diff(rz);

% detect rising edges
rise_edge_idx = find(rz_diff == 1);

% detect falling edges
fall_edge_idx = find(rz_diff == -1);


% There are four possible conditions for running pulses:
% (1) all pulses have the rise edge and fall edge ==>   ____|----|______|----|______|----|_______
% (2) the first pulse does not have the rise edge ==>   ----|______|----|______|----|_______
% (3) the last pulse does not have the fall edge  ==> ____|----|______|----|______|---------
% (4) the first pulse does not have the rise edge and the last pulse
% does not have the fall edge ==> -------|______|----|______|--------

% We treat the conditions in a way as they are condition (4).
% That is we ignore the first and the last pulse of each run.

if length(rise_edge_idx) == length(fall_edge_idx)
    if rise_edge_idx(1) > fall_edge_idx(1) % condition (4)
        rise_edge_idx = [0;rise_edge_idx];
        fall_edge_idx = [fall_edge_idx;0];
    end
else
    if rise_edge_idx(1) > fall_edge_idx(1)
        rise_edge_idx = [0;rise_edge_idx];
    else
        fall_edge_idx = [fall_edge_idx;0];
    end
end

numBlocks = length(rise_edge_idx);
onset = rise_edge_idx + 1;
offset = fall_edge_idx;

for j = 2:numBlocks-1 % the first and the last running blocks are ignored
    idx = idx + 1;
    
    startpoint = offset(1) + 1;
    endpoint = onset(end) - 1;
    
    runningBlocks(idx).session = mouse;
    runningBlocks(idx).mouse_name = mousename;
    runningBlocks(idx).day = day;
    runningBlocks(idx).run = run;
    
    runningBlocks(idx).rz = rz;
    
    runningBlocks(idx).onset = onset(j);
    runningBlocks(idx).offset = offset(j);
    runningBlocks(idx).duration = offset(j) - onset(j);
    runningBlocks(idx).time2previousrun = onset(j) - offset(j-1);
    runningBlocks(idx).time2nextrun = onset(j+1) - offset(j);
    
end

obj.runningBlocks = runningBlocks;

end

