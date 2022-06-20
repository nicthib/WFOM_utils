% Simplified behavior pipeline script. get_behavior() is run for each run in runnames (cell string vector)
function B = behavior_pipeline(runnames)
B = [];
h = waitbar(0,'Getting behavior...');
for n = 1:numel(runnames)
    runname = strrep(runnames{n},'_H.mat','');
    B.(runname) = get_behavior(runnames{n});
    waitbar(n/numel(runnames),h);
end
close(h)
