%%
clear
param_vars = {'sessions_for_analysis','runs_for_analysis'};
DLCdir = '/local_mount/space/dingus/1/RS_analysis/DLC/final';
H_dir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final';
load(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_jan2021','preprocessing_params.mat'),param_vars{:})
addpath(genpath('/local_mount/space/dingus/1/RS_analysis/code'))

for j = 1:numel(runs_for_analysis)
    for i = 1:numel(sessions_for_analysis)
        try
        % session/run info
        file_id = [sessions_for_analysis{i} '_' runs_for_analysis{j}];
        filename = fullfile(H_dir,[file_id '_H.mat']);
        load(filename,'m')
        B = extract_DLC_vars(DLCdir,file_id);
        try m.whisk = sepblockfun(B.whisker2.whisk(1:11980*3),[1 3],@mean); catch; m.whisk = 0; end
        try 
            pB = cleanupDLCvar(B.pupil.Bottom); 
            pT = cleanupDLCvar(B.pupil.Top);  
            m.pupil=pB(:,2)-pT(:,2); 
        catch
            m.pupil = 0; 
        end
        
        % save
        save(filename,'m','-append')
        catch
        end
    end
end
