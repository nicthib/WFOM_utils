clear
load(fullfile('/local_mount/space/dingus/1/RS_analysis/preprocessing','params_92.mat'))
addpath(genpath('/local_mount/space/dingus/1/RS_analysis/code'))
for j = 1:numel(opts.runs_all)
    for i = 1:numel(opts.sessions_all)
        try
        % session/run info
        file_id = [opts.sessions_all{i} '_' opts.runs_all{j}];
        load_dir = findmousefolder(opts.sessions_all{i},opts.runs_all{j},'1');
        mouse_id = opts.sessions_all{i}(1:5);
        m = m_defaults.preprocessing;
        m.IDX = roundIDX(imwarp(m_source.modalIDX,m_all.([mouse_id '_2']).TF,...
            'OutputView',imref2d([m_source.sz m_source.sz]))); 
        m.refim = m_all.([mouse_id '_2']).refim;
        m.baselinewidth = 300; % This is to avoid right after/before running epochs
        m.loadpct = [0 .1]; % For debugging
        m.dofilter = 0; % For debugging
        [m,data,H,Hstd] = LoadData_RS(load_dir,m);
        
        % behavior, summary figure
        B = get_behavior(file_id);
        try m.whisk = B.whisk; catch; m.whisk = 0; end
        try m.pupil = B.pupil; catch; m.pupil = 0; end
        RS_summary_fig(H,m);
        savefig(fullfile(dirs.H,'summary_figs',file_id));
        saveas(gcf,fullfile(dirs.H,'summary_figs',[file_id '.png'])); close all
        
        % remove unused fields to save space
        H = rmfield(H,{'lime','green','red','lime_prefilt','green_prefilt','red_prefilt'});
        try m = rmfield(m,{'aux','DAQ','bkg'}); catch; end % Add a reference about DAQ, aux
        
        % save
        save(fullfile(dirs.H,[file_id '_H.mat']),'H','m')
        catch
        end
    end
end

