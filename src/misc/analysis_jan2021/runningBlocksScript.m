% read data
clear;clc;

% load band-pass filter
addpath /home/bahar/Documents/Auxiliary_Functions
load('/home/bahar/Documents/Filters/lowpassFilter_maximallyFlat.mat')

% latest data files are stored here: 
% addpath /local_mount/space/dingus/1/RS_analysis/analysis_dec2020/H

% for now, I use the old data
addpath /local_mount/space/dingus/1/RS_analysis/H_new

loaddataDIR = '/local_mount/space/dingus/1/RS_analysis/H_new/';

% auxiliary codes
addpath /local_mount/space/dingus/1/RS_analysis/code/utils

mousename = 'cm128';
runs = dir(strcat(loaddataDIR,mousename,'*.mat'));

% remove bad runs
switch mousename
    case 'cm124'
        % remove cm124_1_runI (it's all NaN) and cm124_1_runF (rz is zero)
        runs([4 6]) = [];
    case 'cm125'
        % remove cm125_2_runE (only one principal component is enough to explain
        % 99% variance in the data) and cm125_4_runF (rz is zero)
        runs([10 29]) = [];
    case 'cm126'
        % remove cm126_2_runB, cm126_2_runC, cm126_2_runD
        runs(7:9) = [];
    case 'cm127'
        % remove cm127_1_runH_H.mat (rotf does not look right), cm127_2_runH (rz is zero)
        runs([6 14]) = [];
end

% We use the variable k to follow the number of valid running blocks across
% all days and runs for each mouse. The value of k gets upadated after each
% run.
k = 1;

% We use the variable y to count the number of valid running blocks within
% each run. The value of y gets updated after computing the number of 
% valid running blocks for each run.
y = 0;

% load one day, one run at a time
for n = 1:length(runs)
    
    mouse = runs(n).name;
    mousename = mouse(1:5);
    day = mouse(7);
    run = mouse(9:12);
    
    % load data
    load(mouse)
    
    dims = 1; % chose which number of ROIs you want (1 = 29, 2 = 49, 3 = 83);
    numroi = size(H.jrgeco{dims},1); % number of ROIs
    tnum = size(H.jrgeco{dims},2); % number of time points
    fr = m.framerate; % frame rate of data (x3 for 3 LEDs)
    tt = [0:tnum-1]*3/fr;
    ss = size(H.jrgeco{dims});
    
    % ~6.5 Hz BPF
    neuralsig = H.jrgeco{dims};%(:,1:6000);
    filt = filter(Hd,neuralsig,2);
    
    rz = getrunningpulses(rotf,0.01,20,100);
    
    % detect the running blocks (pulses)
    rz_diff = diff(rz);
    % detect rising edges
    rise_edge_idx = find(rz_diff == 1);
    % detect falling edges
    fall_edge_idx = find(rz_diff == -1);
    
    % examine if the pulses are valid (the criteria is that each valid pulse should have 10 seconds rest befor its onset time and after its offset time)
    T = 10 * 20; % 10 seconds time interval before and after running, which is equivalent to 200 frames
   
    
    % remove the first element of fall_edge idx if it is smaller than the
    % first element of rise_edge idx; this means that a pulse was detected
    % at the begining of running without any rising edge
    if fall_edge_idx(1) < rise_edge_idx(1)
        if length(rise_edge_idx) < length(fall_edge_idx)
            fall_edge_idx = fall_edge_idx(2:end);
        elseif length(rise_edge_idx) == length(fall_edge_idx)
            fall_edge_idx = fall_edge_idx(2:end);
            rise_edge_idx = rise_edge_idx(1:end-1);
        end
    end
    
    % Here, we check if the rise_edge_idx and fall_edge_idx have the same
    % length. If not, we make them have the same length. Sometimes, it's
    % possible that there is no fall edge for a rise edge. This may happen
    % at the end of running. We remove the last index of rise_edge_idx
    if length(rise_edge_idx) > length(fall_edge_idx)
        rise_edge_idx = rise_edge_idx(1:length(fall_edge_idx));
    end
    
    % remove pulses that occured within the first T frames since they do
    % not have T time interval prior to their onset time
    mytemp = rise_edge_idx;
    rise_edge_idx = rise_edge_idx(rise_edge_idx >= T);
    fall_edge_idx = fall_edge_idx(mytemp >= T);
    
    
    num_running_blocks = length(rise_edge_idx);
    
    
    % check if pulses have T rest time interval prior to their onset times and after
    % their offset time
    startpoint = rise_edge_idx - T + 1;
    endpoint = fall_edge_idx + T;
    
    % save the starting points of the valid blocks in T1; save the end points
    % of the valid blocks in T2; save the valid onset times in onset; save the
    % valied offset times in offset
    num_valid_running_blocks = 0;
    for i = 1:num_running_blocks
        % limit the value of endpoint(i) if it exceeds the number of array element
        if endpoint(i) > tnum
            continue
        end
        if (all(~(rz(startpoint(i):rise_edge_idx(i)))) && all(~(rz(fall_edge_idx(i)+1:endpoint(i)))))
            num_valid_running_blocks = num_valid_running_blocks + 1;
            T1(num_valid_running_blocks) = startpoint(i);
            onset(num_valid_running_blocks) = rise_edge_idx(i) + 1;
            offset(num_valid_running_blocks) = fall_edge_idx(i);
            T2(num_valid_running_blocks) = endpoint(i);
        end
    end
    % update the value of y
    y = y + num_valid_running_blocks;
    
    % save the valid blocks along their properties mouse name, day, run, T1,
    % onset time, offset time, T2, neural signal block, running block, pupil
    % block, whisking block
    
    % to go through the indeces of T1 and T2
    j = 0;
    for i = k:y
        j = j + 1;
        runningBlocks(i).onset = onset(j);
        runningBlocks(i).offset = offset(j);
        runningBlocks(i).duration = offset(j) - onset(j);
        runningBlocks(i).neural = filt(:,T1(j):T2(j))';
        runningBlocks(i).rotf = rotf(T1(j):T2(j));
        runningBlocks(i).rz = rz(T1(j):T2(j));
        runningBlocks(i).pupil = pupil(T1(j):T2(j));
        runningBlocks(i).pz = pz(T1(j):T2(j));
        runningBlocks(i).whisk = whisk(T1(j):T2(j));
        runningBlocks(i).wz = wz(T1(j):T2(j));
        runningBlocks(i).mouse_name = mousename;
        runningBlocks(i).day = day;
        runningBlocks(i).run = run;
        runningBlocks(i).runningBlock = j;
    end
    
    % update the value of k
    k = y + 1;
    
    clearvars -except Hd mousename runs k y n runningBlocks
end

save(mousename, 'runningBlocks')




