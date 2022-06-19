%% this function lines up running epochs 
%runnames: list of runnames you want to lineup
%allrunnames: list of total runs that compares directly to st 
%st: states
%frames: number of frames before running
%xlimit: total number of frames

function [runmat] = order_runs(runnames, allrunnames, st, frames, xlimit)


%% Find which indices correspond to st from runnames 
runindex = []
for i = 1:length(runnames)
runindex(i) = find(strcmp(allrunnames, runnames{i})); %% find indices of runnames for st
end

%% Concatonate running and states into 2 large matrices Rot_tot(running) and St_tot(states)
load(runnames{1})
Rot_tot = rotf(1:end-1);
St_tot = st.w2{runindex(1)}';

for i = 2:length(runnames)
    load(runnames{i})
    Rot_tot = horzcat(Rot_tot, rotf(1:end-1));
    St_tot = horzcat(St_tot, st.w2{runindex(i)}');
end
clearvars -EXCEPT runnames St_tot Rot_tot allrunnames st runindex frames xlimit
%% Thresholding: finds where running is above a threshold for a given amount of time
    runmat = zeros(2,xlimit);
    runmat2 = zeros(2,xlimit); %this is just for sorting i'll fix this later

    IsRun = 0; % boolean is running
    r = 0; % running count
    s = 0; % stopping count
    l = 1;
    th = 0.5; % threshold for running 


    for i = frames:length(Rot_tot)-frames % go through all the running data

        if Rot_tot(i) >= th && Rot_tot(i-1) < th && IsRun == 0  %If running is above threshold AND the previous is below AND NOT RUNNING start run epoch
            IsRun = 1;
            ind = i;
            s = 0; % not sure if need this

        elseif IsRun ==1 && abs(Rot_tot(i))>=th %if in running epoch count how long
            r = r +1;
            %s = 0; % start over count
        end

        if abs(Rot_tot(i)) < th && IsRun == 1 %If running is below threshold and IS RUNNING stop count = +1
            s = s +1;
        end
        if s >= 100 && IsRun == 1 && r>20% if stop count is above 20 frames then stop running epoch
            IsRun = 0;
            temp = (xlimit - length(1:r+frames)); % filler 
            runmat(l,:) = St_tot(ind-frames +1:ind+r+temp);
            runmat2(l,1:length(ind-frames +1:ind+r)) = St_tot(ind-frames +1:ind+r);
            l = l+1 ;%index for runmat
            s = 0;
            r = 0;
        end

    end
    %% Sort those running epochs by duration of run
    x = zeros(length(runmat2(:,1)),1);

    for i = 1:length(runmat2(:,1))
        length(find(runmat2(i,:) ~= 0))
        x(i) = length(find(runmat2(i,:) ~= 0));
    end
    [B,I] = sort(x);
    
    runmat = runmat(I,:)

  %figure; imagesc(runmat(I,:))
end
