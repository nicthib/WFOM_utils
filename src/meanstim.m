function stim_data = meanstim(m,data)
% Averages stim responses to a single dataset. Note: m.stimtimes must contain the stimulus times to work. From 2019 RS cohort.
fn = fieldnames(data);

if isfield(m,'stimframes')
    for i = 1:numel(m.stimframes)
        for j = 1:numel(fn)
            if i == 1
                stim_data.(fn{j}) = [];
                stim_data.(fn{j}) = data.(fn{j})(:,:,[-99:200]+m.stimframes(i));
            else
                stim_data.(fn{j}) = stim_data.(fn{j}) + data.(fn{j})(:,:,[-99:200]+m.stimframes(i));
            end
            if i == numel(m.stimframes)
                stim_data.(fn{j}) = stim_data.(fn{j})/numel(m.stimframes);
            end
        end
    end
else
    rz = getbehavioralbouts(smooth(m.rotf,20),.5,20,20);
    rzd = [0; diff(rz)];
    onset_frames = find(rzd == 1);
    offset_frames = find(rzd == -1);
    if onset_frames(1) < 100
        onset_frames(1) = [];
        offset_frames(1) = [];
    end
    n = 0; idx_full = [];
    for i = 1:numel(onset_frames)
        idx = [[-99:20]+onset_frames(i) [0:179]+offset_frames(i)];
        if ~any(m.rotf(idx(150:end))>.5)
            idx_full = [idx_full idx];
            plot(m.rotf(idx)); hold on
            n = n + 1;
            for j = 1:numel(fn)
                if i == 1
                    stim_data.(fn{j}) = [];
                    stim_data.(fn{j}) = data.(fn{j})(:,:,idx);
                else
                    stim_data.(fn{j}) = stim_data.(fn{j}) + data.(fn{j})(:,:,idx);
                end
            end
            if i == numel(onset_frames)
                stim_data.(fn{j}) = stim_data.(fn{j})/n;
            end
        end
    end
end




