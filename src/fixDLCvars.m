% Move pupil variables to "pupil" field
clear
load('pupil_new.mat')
vars = whos;
for i = 1:numel(vars)
    pupil = struct(...
        'Mouse_Run',eval(vars(i).name).Mouse_Run,...
        'Top',eval(vars(i).name).Top,...
        'TopR',eval(vars(i).name).TopR,...
        'Right',eval(vars(i).name).Right,...
        'BottomR',eval(vars(i).name).BottomR,...
        'Bottom',eval(vars(i).name).Bottom,...
        'BottomL',eval(vars(i).name).BottomL,...
        'Left',eval(vars(i).name).Left,...
        'TopL',eval(vars(i).name).TopL,...
        'Eyelid',eval(vars(i).name).Eyelid);
    eval([vars(i).name ' = [];'])
    eval([vars(i).name '.pupil = pupil;'])
    save('pupil_new2.mat',vars(i).name,'-append')
end

%% Combining pupil into H_n files
clear
load(fullfile('/local_mount/space/dingus/1/RS_analysis/DLC','pupil_new2.mat'))
vars = whos;
saveroot = '/local_mount/space/dingus/1/RS_analysis/H_b_new';
for i = 1:numel(vars)
    pupil = eval(vars(i).name).pupil;
    savefile = dir(fullfile(saveroot,['**/' vars(i).name '_H.mat']));
    if ~isempty(savefile)
        save(fullfile(savefile.folder,savefile.name),'pupil','-append')
    end
    i/numel(vars)
end

%% Combining whisk into H_n files
clear
load(fullfile('/local_mount/space/dingus/1/RS_analysis/DLC','whisk.mat'))
clear whisk
vars = whos;
saveroot = '/local_mount/space/dingus/1/RS_analysis/H_b_new';
for i = 1:numel(vars)
    whisk = eval(vars(i).name).whisk;
    savefile = dir(fullfile(saveroot,['**/' vars(i).name '_H.mat']));
    if ~isempty(savefile)
        save(fullfile(savefile.folder,savefile.name),'whisk','-append')
    end
end
