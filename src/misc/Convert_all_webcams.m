% Batch batch convert...
mice = {'cm124','cm125','cm126','cm127','cm128'};
runs = 'BCDEFGHIJ';
savedir = '/local_mount/space/dingus/1/RS_analysis/cm_analysis/train_webcam';
for a = 1:numel(mice)
    for b = 1:8
        mouse = [mice{a} '_' mat2str(b)];
        for c = 1:numel(runs)
            mdir = findmousefolder(mouse);
            mdirfull = fullfile(mdir,'webcam',['run' runs(c)]);
            vf = 1:numel(dir(mdirfull))/2-1;
            try
                BatchLoadFLIR(mouse,['run' runs(c)],vf,60,savedir,0)
            catch
                disp(['failed run ' runs(c) ' for mouse ' mice{a} '_' num2str(b)])
            end
        end
    end
end