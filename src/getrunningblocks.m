function runningBlocks = getrunningblocks(r,T,runname)
rz = getbehavioralbouts(r.rotf,1,1,1); % .5 for non-smoothed, .1 for smoothed?
b_bool = getbehavioralbouts(r.rotf,1,10,10);
L = numel(rz);

% detect the running blocks (pulses)
b_diff = diff(b_bool);
% detect rising edges
rise_edge_idx = find(b_diff == 1);
% detect falling edges
fall_edge_idx = find(b_diff == -1);

[ttn,ttp] = behav_time(rz);

% There are four possible conditions for running pulses:
% (1) all pulses have the rise edge and fall edge ==>   ____|----|______|----|______|----|_______
% (2) the first pulse does not have the rise edge ==>   ----|______|----|______|----|_______
% (3) the last pulse does not have the fall edge  ==> ____|----|______|----|______|---------
% (4) the first pulse does not have the rise edge and the last pulse
% does not have the fall edge ==> -------|______|----|______|--------

% We treat the conditions in a way as they are condition (4).
% That is we ignore the first and the last pulse of each run.

% number of blocks to consider
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
idx = 1;
for j = 2:numBlocks-1 % the first and the last running blocks are ignored
    if onset(j) - T(1) < 0 % ignore the running blocks with the onset less than T frames
        continue;
    else
        startpoint = onset(j) - T(1);
    end
    
    if onset(j+1) + T(2) > L % ignore the running blocks with offsets less than T frames before the end of the session
        continue;
    else
        endpoint = offset(j) + T(2); % Each 'running block' is now from onset to onset.
    end
    
    runningBlocks(idx).onset = onset(j); % onset(j) + 1 is the actual run
    runningBlocks(idx).offset = offset(j);
    runningBlocks(idx).offsetduration = onset(j+1) - offset(j);
    runningBlocks(idx).onsetduration = offset(j) - onset(j);
    runningBlocks(idx).time2previous = onset(j) - offset(j-1);
    runningBlocks(idx).time2next = onset(j+1) - offset(j);
    
    runningBlocks(idx).d = r.d(startpoint:endpoint,:);
    runningBlocks(idx).rotf = r.rotf(startpoint:endpoint);
    
    runningBlocks(idx).pupil = hampel(r.pupil(startpoint:endpoint),50);
    
    runningBlocks(idx).whisk = r.whisk(startpoint:endpoint);
    
    runningBlocks(idx).p_score = r.p_score;
    
    runningBlocks(idx).ttp = ttp(startpoint:endpoint);
    runningBlocks(idx).ttn = ttn(startpoint:endpoint);
    
    runningBlocks(idx).runname = runname;
    idx = idx + 1;
end
%disp(['found ' mat2str(numBlocks) ' bouts'])
