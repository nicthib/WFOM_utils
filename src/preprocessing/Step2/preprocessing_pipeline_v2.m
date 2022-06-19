clear
load(fullfile('/local_mount/space/dingus/1/RS_analysis/preprocessing','params_92.mat'))
addpath(genpath('/local_mount/space/dingus/1/RS_analysis/code'))
dirs.savepath = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_92_RAW';
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
        m.dofilter = 0;
        m.outputs = 'rgl';
        [m,data,H] = LoadData_RS_RAW(load_dir,m);
        m.checkim = data.lime(:,:,100);
        %imshowpair(m.checkim,logical(imgradient_nic(m.IDX))); axis image; axis off;
        
        % save
        save(fullfile(dirs.savepath,[file_id '_H.mat']),'H','m')
        catch
        end
    end
end

