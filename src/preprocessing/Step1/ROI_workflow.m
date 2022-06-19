%% ROI map script
% This script seeks to create ROI maps using a semi-supervised clustering
% approach. Each map is adjusted on a subject level using a data-driven
% approach, and is composed of bilaterally matching components. 

% There are n steps to this process:
% 1. Source map creation
% 2. Registration
% 3. Secondary clustering and creation of modal cluster assignment

%% 1a. Source map creation
% Here, we load a single run to create a source map, which will be used as
% an initial condition for generating all ROI maps. We use cm125_2_run_D
% for this cohort.
%%
m = m_defaults.preprocessing; % Default options
m.outputs = 'Pg'; m.PCAcomps = 500; % PCA for smaller file size
[m,data] = LoadData_v2(findmousefolder(opts.source_session,opts.source_run,1),m);
m.refim = data.green(:,:,1); % refim here refers to a single green image used to register ROI maps.
data = rmfield(data,'green'); % We only needed green for the refim
m = rmfield(m,{'aux','DAQ'}); % Saves a bit more space
save(fullfile('/local_mount/space/dingus/1/RS_analysis/preprocessing/Step1','source_data.mat'),'data','m','-v7.3')
% We can now load this later when needed.

%% 1b. Reconstructing source map dataset
clear
load('/local_mount/space/dingus/1/RS_analysis/preprocessing/Step1/source_data.mat')
load('params_checkpoint.mat')
data = recon_SVD(data.P,m); % Load the data we previously saved and reconstruct it from the SVD components
jrgeco = data.jrgeco; clear data; % clear up space
% Now you have jrgeco and m for the next step

%% 1c. Creating map using GUI
% Run GUI here and select components. The goal is to create a map with 
% an equal number of bilaterally matching components. 
m_source = IDX_GUI(m,jrgeco);

%% Save checkpoint file
save(fullfile(dirs.param_file,'params_checkpoint.mat'),'m_source','-append')
% m_source.IDX_final is now the starting point for all ROI maps.

%% 2a. Getting affine transforms for each mouse
% Here, we load a single image from each mouse, and run another GUI that
% allows us to align the source map to the given dataset. This GUI returns
% the affine transform needed to translate the source map to a given
% dataset.
for i = 1:numel(opts.sessions_registration)
    fn = opts.sessions_registration{i}; % Session name
    m = m_defaults.preprocessing;
    m.loadpct = [0 .01]; m.outputs = 'g'; % We only need one image, so we just load 1%
    [~,data] = LoadData_v2(findmousefolder(opts.sessions_registration{i}...
        ,opts.runs_registration{i},1),m); % Load data here
    m_all.(opts.sessions_registration{i}).refim = data.green(:,:,1); % The first frame will be used as a refernce image
    m_all.(fn).TF = TF_GUI(m_source.IDX_LR_refined,m_all.(fn).refim); % Run GUI here for each session. Close button will return the TF
    figure; imagesc(m_all.(opts.sessions_registration{i}).refim); axis image; axis off
    [x,y] = getpts; close all; % Click two points that form the midline
    m_all.(fn).midline = [x y]; % Midline will be used later to split L and R hemispheres
end

% We now have transformation matrices to map the source map to each
% session.
save(fullfile(dirs.param_file,'params_checkpoint.mat'),'m_all','-append')

%% 2b. Show results (optional)
close all
f1 = figure('Position',[0 0 1000 1000]);
f2 = figure;
for n = 1:numel(opts.sessions_registration)
    fn = opts.sessions_registration{n};
    gtmp = m_all.(fn).refim;
    im_TF = imwarp(m_source.IDX_LR_refined,m_all.(fn).TF,'OutputView',imref2d([m_source.sz m_source.sz]));
    figure(f1);
    imagesc(gtmp); colormap gray; axis image;
    im = getframe; im = imresize(im.cdata,[m_source.sz m_source.sz]);
    figure(f2); subplot(2,3,n)
    imagesc(imoverlay(im,imerode(logical(imgradient_nic(im_TF)),ones(2)),'w'))
    axis image; axis off
end
close(f1)

%% 3a. Secondary clustering
% In this section, we load datasets from each mouse, and calculate kmeans
% maps from 10 second epochs, using the m_source map as an initial
% condition. This is similar to using replicates in traditional clustering,
% but we want to avoid using random initialization due to the cluster label
% inconsistencies. Once we have calculated ROI maps from each epoch, we
% register these back to the source map, and calculate a final ROI map
% using the most common cluster assignments for each pixel in the FOV. We
% can then use this as a "master map" for all TC extraction.

%%
% modalIDX is the most common cluster assignments across numerous epochs
load('params_checkpoint.mat')
% m_source.IDX_final is the m_source map. First, we "unilateralize" it.
m_source.IDX_final = m_source.IDX_LR_refined;
nU = max(m_source.IDX_final(:))/2; % # of unilateral components (should be 46)
idx = find(m_source.IDX_final > nU); % Find components that are from teh right side
m_source.IDX_final(idx) = m_source.IDX_final(idx)-nU; % Subtract by nU to make a unilateral map with 46 components
m_source.IDX_finalL = m_source.IDX_final.*m_source.L; % Create L and R IDX maps
m_source.IDX_finalR = m_source.IDX_final.*m_source.R;
%IDX1_n = IDXtodigraph(m_source.IDX_final); % Old approach (by comparing
% neighbors instead of modal approach)
allIDXs = []; % allIDXs is each kmeans output across all subjects.
for i = 1:numel(opts.sessions_registration)
    fn = opts.sessions_registration{i};
    m = m_defaults.preprocessing;
    fns = fieldnames(m_all.(fn)); % Combine m_all with default m for loading data
    for f = 1:numel(fns); m.(fns{f}) = m_all.(fn).(fns{f}); end
    loadpath = findmousefolder(opts.sessions_registration{i},opts.runs_registration{i},'1');
    [m,data] = LoadData_v2(loadpath,m);
    rotf = smooth(abs(m.rotf),200);
    m.score = []; m.epochs = [];
    % Now, we get maps for each side of the brain (since this is how we got
    % the source map in the first place);
    IDX_TFL = roundIDX(imwarp(m_source.IDX_finalL,m_all.(fn).TF,...
        'OutputView',imref2d([m_source.sz m_source.sz])));
    IDX_TFR = roundIDX(imwarp(m_source.IDX_finalR,m_all.(fn).TF,...
        'OutputView',imref2d([m_source.sz m_source.sz])));

    IDX_TFL(IDX_TFL==0) = NaN; 
    IDX_TFR(IDX_TFR==0) = NaN;
    BWL = double(IDX_TFL>0); 
    BWL(BWL==0) = NaN;
    BWR = double(IDX_TFR>0); 
    BWR(BWR==0) = NaN;
    % Choose epochs with no running
    epoch_length = 200; % 10 second epochs
    for e = 1:epoch_length:numel(rotf) - epoch_length
        if max(rotf(e:e+epoch_length)) < .5
           m.epochs(end+1) = e;
        end
    end
    disp(['found ' mat2str(numel(m.epochs)) ' epochs'])
    
    % Loop through each epoch, extract H, run kmeans with these initial
    % conditions, and then append the result to the rest.
    for e = 1:numel(m.epochs)
        kmeans_idx = m.epochs(e):m.epochs(e)+epoch_length;
        jrgecoL = data.jrgeco(:,:,kmeans_idx).*repmat(BWL,[1 1 numel(kmeans_idx)]);
        jrgecoR = data.jrgeco(:,:,kmeans_idx).*repmat(BWR,[1 1 numel(kmeans_idx)]);
        HL = getHfromKmeans(jrgecoL,IDX_TFL,1);
        HR = getHfromKmeans(jrgecoR,IDX_TFR,1);
        
        IDXL = reshape(kmeans(reshape(jrgecoL,[256^2 size(jrgecoL,3)]),size(HL,1),'Distance','correlation','Start',HL),[m_source.sz m_source.sz]);
        IDXR = reshape(kmeans(reshape(jrgecoR,[256^2 size(jrgecoR,3)]),size(HR,1),'Distance','correlation','Start',HR),[m_source.sz m_source.sz]);
        IDXLR = nansum(cat(3,IDXL,IDXR),3);
        m.IDXs(:,:,e) = IDXLR;
        TFi = m_all.(fn).TF; T = pinv(TFi.T); T(:,3) = [0 0 1]; TFi.T = T;
        allIDXs(:,:,end+1) = roundIDX(imwarp(m.IDXs(:,:,e),TFi,'OutputView',imref2d([m_source.sz m_source.sz])));
        %IDX2_n = IDXtodigraph(m.IDXs(:,:,e)); % Old code
        %try m.score(e) = graph_similarity(IDX1_n,IDX2_n); catch; end % Old code
    end
    m.modalIDX = mode(m.IDXs,3);
    try m = rmfield(m,{'DAQ','aux'}); catch; end
    m_all.(fn) = m;
end

%% Subject level modal ROI maps
close all;
figure('Position',[0 0 1000 750])
ha = tight_subplot(2,3,[0 0],[.1 .1],[0 0]);
for i = 1:numel(opts.sessions_registration)
    fn = opts.sessions_registration{i};
    [L,R] = midline(m_all.(fn).midline);
    IDX_U = m_all.(fn).modalIDX;
    m_all.(fn).IDX_final = IDX_U+(double(R.*(IDX_U>0)).*max(IDX_U(:)));
    axes(ha(i))
    n = max(m_all.(fn).IDX_final(:));
    showIDX_labels_subregions(m_all.(fn).IDX_final,...
        networks,...
        [lines(numel(networks))],1,ha(i));
    title(strrep(opts.sessions_registration{i},'_',' '))
end

%% Overall Modal ROI map
IDXtmp = cleanupIDX(mode(allIDXs,3),50);
IDXtmp(m_source.BW==0) = NaN;
modalIDX = nnIDX(IDXtmp);
networks = max(modalIDX(:));
IDXscore = mean(allIDXs==mode(allIDXs,3),3);
IDXscore(IDXscore==1) = NaN;
figure; 
subplot(121)
imagesc(IDXscore); axis image; axis off; caxis([.5 1])
colormap(jet)
im = getframe(gca); im = imresize(im.cdata,size(modalIDX));
imagesc(imoverlay(im,imerode(logical(imgradient_nic(modalIDX)),ones(2)),'w'))
axis image; axis off;
colorbar
title('Modal assignment score')
caxis([.5 1])

subplot(122)
n = max(m_source.IDX_final(:))*2;
showIDX_labels_subregions(modalIDX,networks,lines(numel(networks)),1,gca);
title('Modal ROI map')

%% Final cleanup, if needed
close all;
showIDX_labels_subregions(modalIDX,networks,lines(numel(networks)),0,gca);
[x,y] = getpts; x = round(x); y = round(y);
comp = modalIDX(y,x);
cc = bwconncomp(modalIDX==comp);
for i = 1:numel(cc.PixelIdxList)
   ROI = zeros(size(modalIDX));
   ROI(cc.PixelIdxList{i}) = 1;
   if ROI(y,x) == 1
       roiofinterest = cc.PixelIdxList{i};
   end
end
newcomp = inputdlg('New value:');
modalIDX(roiofinterest) = str2num(newcomp{1});
close all;
showIDX_labels_subregions(modalIDX,networks,lines(numel(networks)),0,gca);

%% If good, save
m_source.modalIDX = modalIDX - m_source.BW_R*max(modalIDX(:));
save(fullfile(dirs.param_file,'params_checkpoint.mat'),'m_all','m_source','-append')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OLD CODE FOR REGISTRATION (SUCKS)
% clear; load('params_checkpoint.mat')
% m_source.cp = [128,50;128,206]; % Initial cps
% m_source.refim = uint8(m_source.refim*255/max(m_source.refim(:)));
% for i = 1:numel(opts.sessions_registration)
%     % Need to update
%     m = m_defaults.preprocessing;
%     m.loadpct = [0 .01]; m.outputs = 'g';
%     [~,data] = LoadData_v2(findmousefolder(opts.sessions_registration{i}...
%         ,opts.runs_registration{i},1),m);
%     data.green = data.green(:,:,1);
%     
%     % ms is the mouse of interest
%     m_all.(opts.sessions_registration{i}).refim = ...
%         uint8(data.green*255/max(data.green(:)));
%     m_all.(opts.sessions_registration{i}).cp = m_source.cp;
%     
%     % Select registration points
%     % Here, we show the points on source map, and you are expected
%     % to adjust the points on the right image to match. We use this to
%     % get a transformation from one mouse to the source map.
%     [m_source.cp,m_all.(opts.sessions_registration{i}).cp] = cpselect(...
%     m_source.refim,...
%         m_all.(opts.sessions_registration{i}).refim,...
%         m_source.cp,m_source.cp,'Wait',true);
%   
%     % Register masks
%     m_all.(opts.sessions_registration{i}).TF = ...
%         fitgeotrans(m_source.cp,m_all.(opts.sessions_registration{i}).cp...
%         ,'NonreflectiveSimilarity');
% end
