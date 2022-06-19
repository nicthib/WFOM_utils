% Step 1: BW masks for getting TF matrices
% Default m struct and m_all for later
load(fullfile('/local_mount/space/dingus/1/RS_analysis/Draft','params.mat'),'m_default','m_all')
ROImap = m_all.cm125.IDX_final; % The ROI map we will eventualy transform for other mice
clear m_all
% This list can be added to as more mice are added
session_names = {'cm124_2','cm125_2','cm126_2','cm127_2','cm128_2','cm193_1'};
run_names = {'runD','runD','runD','runD','runD','runA'};
for s = 1:numel(session_names)
    session_name = session_names{s}; run_name = run_names{s};
    % Where we are loading data from
    load_dir = findmousefolder(session_name,run_name,'1');
    % First just load aux data to select a kmeans epoch
    m = m_default; m.outputs = 'gn'; m.noload = 1;
    % cm193_1 had LED switching issues, and I also don't know the bkg
    if strcmp(session_name,'cm193_1')
        m.LED_offset = 2; m.bkgsub = [];
    end
    [m,~] = LoadData_v21(load_dir,m);
    % Pick kmeans epoch (no running)
    close all; plot(m.rotf); [x,~] = ginput(2); x = round(x); close all
    % Now load that epoch
    m.loadpct = x/numel(m.rotf); m.noload = 0;
    m.PCAcomps = 0; m.BW = ones(256,256); m.baseline = 1:100; % Reset baseline
    [m,data] = LoadData_v21(load_dir,m);
    jrgeco = reshape(data.jrgeco,[256^2 size(data.jrgeco,3)]);
    
    % Now, make BW mask
    IDX = reshape(kmeans(jrgeco,100,'distance','correlation'),[256 256]);
    IDXtmp = zeros(256,256);
    % Clean up by removing objects less than 200 pix and eroding 
    for i = 1:max(IDX(:))
        tmp = IDX==i; tmp = imerode(imfill(bwareaopen(tmp,200),'holes'),ones(3));
        IDXtmp(tmp==1) = 1;
    end
    % Get conv hull of result, erode again
    IDXtmp = bwconvhull(IDXtmp);
    IDXtmp = imerode(IDXtmp,ones(4));
    jrgeco(IDXtmp(:)==0,:) = NaN;
    % Recluster
    IDX = reshape(kmeans(jrgeco,100,'distance','correlation'),[256 256]);
    
    % Final crop. Make sure each one has a consistent shape so registration
    % is successful
    imagesc(IDX); axis image; axis off;
    BWs.([session_name '_' run_name]) = roipoly;
    close all
end
% Only important variable here is BWs. That is needed for next step

%% Assign to IDXstart for next step, show and check
moving = logical(ROImap); % Needs to be put in dedicated variable
fns = fieldnames(BWs);
for i = 1:6
    session_name = session_names{i};
    run_name = run_names{i};
    m = m_default; m.outputs = 'l'; m.loadpct = [0 .001]; m.BW = ones(256,256);
    load_dir = findmousefolder(session_name,run_name,'1');
    if strcmp(session_name,'cm193_1')
        m.LED_offset = 2; m.bkgsub = [];
    end
    [m,data] = LoadData_v21(load_dir,m);
    fixed = BWs.(fns{i});
    subplot(2,6,i)
    imshowpair(moving,fixed)
    [optimizer, metric] = imregconfig('monomodal');
    TF = imregister_nic(double(moving), double(fixed), 'affine', optimizer, metric);
    
    subplot(2,6,i+6)
    IDX_start.(fns{i}(1:7)) = roundIDX(imwarp_fr(ROImap,TF));
    imshowpair(imwarp_fr(moving,TF),fixed)
end
% Check the top row (before) and compare to the bottom (after). There
% should be a clear improvement in overlap.

%% Reclustering
load(fullfile('/local_mount/space/dingus/1/RS_analysis/Draft/params.mat'),'m_default')
session_names = {'cm124_2','cm125_2','cm126_2','cm127_2','cm128_2','cm193_1'}; % cm125_2 runD is the base map, no need to recalculate it.
run_names = {'runD','runD','runD','runC','runD','runA'};
for run_i = 1:5
    file_id = [session_names{run_i} '_' run_names{run_i}];
    load_dir = findmousefolder(session_names{run_i},run_names{run_i},'1');
    m = m_default; m.outputs = 'gn';
    if strcmp(session_names{run_i},'cm193_1')
        m.LED_offset = 2; m.bkgsub = [];
    end
    m.loadpct = [0 1]; m.BW = ones(256,256);
    [m,data] = LoadData_v21(load_dir,m);
    m.aux = []; m.DAQ = [];
    rotf = smooth(abs(m.rotf),200);
    m.score1 = []; m.score2 = [];  m.score3 = []; m.epochs = [];

    % remove epochs with running
    epoch_length = 1200;
    for i = 1:200:10000
        if max(rotf(i:i+1200)) < .5
           m.epochs(end+1) = i;
        end
    end
    m.refim = data.green(:,:,1);

    for i = 1:numel(m.epochs)
        m.kmeans_idx = m.epochs(i):m.epochs(i)+epoch_length;
        jrgeco = data.jrgeco(:,:,m.kmeans_idx);

        IDX1 = IDX_start.(session_names{run_i});
        IDX1(IDX1==0) = NaN;
        BW = double(IDX1>0); BW(BW==0) = NaN;
        BW_R = double(IDX1 <= 92/2); BW_R(BW_R == 0) = NaN;
        BW_L = double(IDX1 >  92/2); BW_L(BW_L == 0) = NaN;
        jrgeco = jrgeco.*repmat(BW,[1 1 size(jrgeco,3)]);
        m.H = getHfromKmeans(jrgeco,IDX1,1);
        m.IDX(:,:,i) = reshape(kmeans(reshape(jrgeco,[256^2 size(jrgeco,3)]),size(m.H,1),'Distance','correlation','Start',m.H),[256 256]);
        IDX2 = m.IDX(:,:,i);

        % Metrics
        % Check overlap
        denom = sum(IDX2(:)>=1);
        numer = sum(abs(IDX1(:)-IDX2(:))==0);
        m.score1(i) = round(numer/denom,2);

        % Check graph similarity
        IDX1_n = IDXtodigraph(IDX1);
        IDX2_n = IDXtodigraph(IDX2);
        m.score2(i) = round(graph_similarity(IDX1_n,IDX2_n),3);

        % Check centroid nearness
        m.score3(i) = round(centroid_nearness(IDX1,IDX2),1);
        % Display scores
        disp(['Score1: ' mat2str(m.score1(i)) ', Score2: ' mat2str(m.score2(i)) ', Score3: ' mat2str(m.score3(i))])
    end
    besti = find(m.score2 == max(m.score2));
    
    % Assign to results
    m.IDX_best = m.IDX(:,:,besti);
    m_all.(m.mouse) = m;
end

%% Reset to 46 component map, then split along midline and group back to 92 component map.
for i = 6
    m = m_all.(session_names{i});
    if i == 2 % skip for base map
       m.IDX_final = IDX_start.cm125_2;
       m_all.(session_names{i}) = m;
       continue
    end
    IDX = m.IDX_best;
    close all
    imagesc(logical(imgradient_nic(IDX))); axis image; axis off
    title('Draw a boundary around the left hemisphere')
    BW_L = roipoly;
    IDX(IDX>46) = IDX(IDX>46) - 46;
    m.IDX_final = IDX.*BW_L+((IDX+46).*~BW_L);
    m_all.(session_names{i}) = m;
end

%% fix cm127
% IDX = m_all.cm127_2.IDX_final;
% cand = [9 35];
% imagesc(ismember(IDX,[cand cand+46])); axis image; axis off
% BW = roipoly;
% IDX(and(BW,ismember(IDX,[cand]))==1) = cand(1);
% IDX(and(BW,ismember(IDX,[cand+46]))==1) = cand(1)+46;
% IDX(and(~BW,ismember(IDX,[cand]))==1) = cand(2);
% IDX(and(~BW,ismember(IDX,[cand+46]))==1) = cand(2)+46;
% m_all.cm127_2.IDX_final = IDX;

%% Save and show
% clearvars -EXCEPT m_all IDX_start
% save('IDX_final.mat')
clear
load('/local_mount/space/dingus/1/RS_analysis/Draft/params.mat','final_ord','m_all')
figure('Position',[0 0 1920 500])
session_names = {'cm124','cm125','cm126','cm127','cm128'};
ha = tight_subplot(1,5,[0 0],[.1 .1],[0 0]);
for i = 1:5
    axes(ha(i))
    m = m_all.(session_names{i});
    IDX = m.IDX_final;
    IDX = remap_IDX(IDX,final_ord);
    showIDX_labels_subregions(IDX,[8 8 7 4 7 12],lines(6),1,gca);
    title(strrep(session_names{i}(1:5),'_',' '))
end
%saveas(gcf,'IDX_final.fig')
%saveas(gcf,'IDX_final.png')

%% Comparing FC maps
clearvars -EXCEPT m_all IDX_start
load(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_jan2021','preprocessing_params.mat'),'m_default')
fns = fieldnames(m_all);
for i = 1:numel(fns)
    m = m_all.(fns{i}).m; m.PCAcomps = 0; m.IDX = m.IDX_final;
    [m,data] = LoadData_v2(m.fulldir,m);
    H_all.(fns{i}) = data.H.jrgeco;
end

%%
fns = {'cm124_2','cm125_2','cm126_2','cm127_2','cm128_2'};
regions = {'M-L','SS-L','SSWA-L','RS-L','V-L','M-R','SS-R','SSWA-R','RS-R','V-R'};
for i = 1:5
    subplot(2,3,i)
    cc(:,i) = reshape(corr(H_all.(fns{i})(:,baselinefromrot(m_all.(fns{i}).rotf,1000,.5))'),[92^2 1]);
    show_state_centroid(reshape(cc(:,i),[92 92]),[9 16 9 5 7],regions);
    title(fns{i}(1:5))
end
subplot(236)
imagesc(triu(corr(cc))-diag(diag(corr(cc))))
axis image
cmap = jet;
cmap(1,:) = 0;
colormap(cmap)
colorbar
caxis([0 1])
xticks(1:5); yticks(1:5)
xtickangle(45)
xticklabels(cellfun(@(x) x(1:5),fns,'UniformOutput',false))
yticklabels(cellfun(@(x) x(1:5),fns,'UniformOutput',false))
