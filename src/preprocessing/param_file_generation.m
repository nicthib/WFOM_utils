%% Step 0: Some directories, filenames, and metadata
% Directories
clear
dirs = [];
dirs.DLC = '/local_mount/space/dingus/1/RS_analysis/DLC/final'; % Where DLC outputs are located
dirs.savepath = '/local_mount/space/dingus/1/RS_analysis/preprocessing';
dirs.H = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_92comps_new';

structfun(@(x) (addpath(genpath(x))),dirs);

% Sessions, runs
opts = [];
opts.sessions_all = {'cm124_1','cm124_2','cm125_1','cm125_2','cm125_3','cm125_4',...
    'cm125_5','cm125_6','cm125_7','cm126_1','cm126_2','cm126_3',...
    'cm126_5','cm126_6','cm127_1','cm127_2','cm127_3','cm128_1',...
    'cm128_2','cm128_3','cm128_4','cm128_5','cm128_6','cm128_7','cm128_8'};
% All (possible) run names
opts.runs_all = {'runB','runC','runD','runE','runF','runG','runH','runI','runJ'};
% One session from each mouse to obtain registration
opts.sessions_registration = {'cm124_2','cm125_2','cm126_2','cm127_2','cm128_2'};
% One run from each above session for registration
opts.runs_registration = {'runD','runD','runD','runD','runD'};
opts.source_session = 'cm125_2';
opts.source_run = 'runD';

% Filter
load('/local_mount/space/dingus/1/RS_analysis/preprocessing/filter.mat')
% Hd = []; Developed by Bahar using the filter GUI (?)

%% Default m's
m_defaults = [];
% For correction analysis
m = makem; m.loadpct = [0 .5]; m.outputs = 'rgl';
m.BW = ones(256,256); m.dsf = 2; m.nrot = 2;
m_defaults.correction_analysis = m;

% preprocessing
m = makem; % Create a premade m to save time
m.outputs = 'rglodn'; % We want all outputs so we can plot some raw TCs on the summary figure
m.loadpct = [0 .9995]; % This is done to truncate the final spool, which contains 1 real frame of data and 4 frames of nothing. Avoids dealing with NaNs.
m.dsf = 2; % We downsample to 256x256 images
% Dr,Dg were estimated from code located in: /local_mount/space/dingus/1/RS_analysis/preprocessing 
% See the file Mean_obj_function.png for the general range of acceptable values.
m.Dr = 1.1; m.Dg = .2;
m.baselinewidth = 300; % 300 frames are initially picked, then the middle 100 are finally used as a baseline. This is to avoid using frames right before/after movement.
m.nrot = 2; % To rotate FOV to proper orientation
m.bkgsub = [];
m.spectra_file = '/local_mount/space/dingus/1/RS_analysis/spectra/spectra_010621_WFOM2.mat'; % Spectra file for hemodynamic conversion
m.BW = ones(256,256); % We do not crop, since we use an IDX map for TC extraction
m.dofilter = 1;
m.Hd = Hd;
m_defaults.preprocessing = m;
clear m i fns n

%% Save checkpoint
save(fullfile(dirs.param_file,'params_checkpoint.mat'))

%% Step 1: ROI maps and transformation matrices
% Run code in the step 1 folder. At the end, you should have 'm_source' and
% 'm_all'

%% Step 2: H ordering and networks (optional)
% Run code in the step 2 folder. At the end, you should have `final_ord` and
% `networks` added to the params file.

%% Step 3 - final runname pointers
allrunnames = getallexps(dirs.H,{'H'})'; % Gets all runs in the folder with the criteria above
runnames_RS = []; runnames_AS = []; % Resting State and "Active State" runs should be split. 
for i = 1:numel(allrunnames)
    load(fullfile(dirs.H,allrunnames{i}),'m')
    if isfield(m,'stimframes')
        runnames_AS{end+1} = allrunnames{i};
    else
        runnames_RS{end+1} = allrunnames{i};
    end
end
runnames_AS = runnames_AS'; runnames_RS = runnames_RS';

save('params_checkpoint.mat','runnames_AS','runnames_RS','-append')
