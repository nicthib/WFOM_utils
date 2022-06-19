function results = get_nls_frames(H,behavior,ww)

% Set default options
f = 1;
for i = 1:ww:size(H,2)-ww
    idx = i:i+ww;
    if ~any(isnan(behavior.pupil(idx)))
        results.n(f,:) = reshape(corr(H(:,idx)'),[92^2 1]);
        results.rotf(f,:) = behavior.rotf(idx);
        results.whisk(f,:) = behavior.whisk(idx);
        results.pupil(f,:) = behavior.pupil(idx);
        f = f+1;
    end
end




