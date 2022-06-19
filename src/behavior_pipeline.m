function B = behavior_pipeline(runnames)
B = [];
h = waitbar(0,'Getting behavior...');
for n = 1:numel(runnames)
    runname = strrep(runnames{n},'_H.mat','');
    B.(runname) = get_behavior(runnames{n});
    waitbar(n/numel(runnames),h);
end
close(h)
%save(fullfile('/local_mount/space/dingus/1/RS_analysis/Draft/dFC','behavior.mat'),'B','-v7.3')
