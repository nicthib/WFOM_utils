clearvars -EXCEPT IDX_final m_all
param_vars = {'m_default','IDX_registered'};
load(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_jan2021','preprocessing_params.mat'),param_vars{:})
session_names = {'cm126_2','cm127_2','cm128_2'}; % cm125_2 runD is the base map, no need to recalculate it.
run_names = {'runD','runD','runD','runD'};
run_i = 3;
file_id = [session_names{run_i} '_' run_names{run_i}];
load_dir = findmousefolder(session_names{run_i},run_names{run_i},'1');
m = m_default; m.outputs = 'gn'; %m.noload = 1;
m.loadpct = [0 .2]; m.PCAcomps = 500; m.BW = ones(256,256); 
[m,data] = LoadData_v2(load_dir,m);

%%
close all; plot(m.rotf); 
[kmeans_idx,~] = ginput(2); kmeans_idx = round(kmeans_idx);
close all
jrgeco = data.jrgeco(:,:,kmeans_idx(1):kmeans_idx(2));
m.kmeans_idx = kmeans_idx;

% Make an exclusion mask if necessary
BW_exclude = zeros(256,256); BW_exclude = jrgeco(:,:,50)==-1;

% Split L and R
IDX1 = IDX_registered.(session_names{run_i});
IDX1(BW_exclude==1) = NaN;
IDX1(IDX1==0) = NaN;
BW_R = double(IDX1 <= 92/2); BW_R(BW_R == 0) = NaN;
BW_L = double(IDX1 >  92/2); BW_L(BW_L == 0) = NaN;
jrgeco_L = jrgeco.*repmat(BW_L,[1 1 size(jrgeco,3)]);
jrgeco_R = jrgeco.*repmat(BW_R,[1 1 size(jrgeco,3)]);
H = getHfromKmeans(jrgeco,IDX1,0);
H_L = H(1:92/2,:); H_R = H(92/2+1:end,:);

IDX_L = kmeans(reshape(jrgeco_L,[256^2 size(jrgeco,3)]),46,'Distance','correlation','Start',H_L+H_R);
IDX_R = kmeans(reshape(jrgeco_R,[256^2 size(jrgeco,3)]),46,'Distance','correlation','Start',H_R+H_L);

IDX2 = reshape(nansum([IDX_L+46 IDX_R]'),[256 256]);
IDX2(IDX2==0) = NaN;

% Check overlap
denom = sum(IDX2(:)>=1);
numer = sum(abs(IDX1(:)-IDX2(:))==0);
score = numer/denom;


% Make cIDX, show
subregions = [5 9 9 7 16];
subplot(1,2,1)
showIDX_labels_subregions(IDX1,subregions) 
subplot(1,2,1)
showIDX_labels_subregions(IDX2,subregions) 

% Assign
IDX_final.(session_names{run_i}) = IDX2;
m_all.(session_names{run_i}) = m;





%% Register basemaps to maps
load(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_jan2021','preprocessing_params.mat'),'IDX')
basemap = IDX2.cm125_2;
IDX2 = rmfield(IDX2,'cm125_2');
IDX_registered.cm125_2 = basemap;
fn = fieldnames(IDX2);
for i = 1:numel(fn)
    [optimizer, metric] = imregconfig('monomodal');
    fixed = IDX2.(fn{i});
    moving = basemap;
    TF = imregister_nic(double(moving >= 1), double(fixed >= 1), 'affine', optimizer, metric);
    IDX_registered.(fn{i}) = roundIDX(imwarp_fr(basemap,TF));
end
save(fullfile('/local_mount/space/dingus/1/RS_analysis/analysis_jan2021','preprocessing_params.mat'),'IDX_registered','-append')

%% Show final
session_names = {'cm124_2','cm125_2','cm126_2','cm127_2','cm128_3'};
subregions = [9 16 9 5 7];
close all; figure('Position',[0 0 1200 1000])
ha = tight_subplot(2,3,[.02 .02],[.1 .1],[.01 .01]);
for i = 1:5
    fn = session_names{i};
    IDXtmp = IDX_erode.(fn);
    axes(ha(i));
    showIDX_labels_subregions(IDXtmp,subregions) 
    title(strrep(fn,'_',' '))
end

%% New order

ord = [10 9 8 7 11 12 13 14 6, ...
40 37 36 39 38 35 31 33 32 42 41 34 43 45 44 46,...
15 16 17 18 19 20 21 22 23,...
3 5 2 4 1,...
27 26 29 25 24 28 30];

%% Flipping L and R

session_names = {'cm124_2','cm125_2','cm126_2','cm127_2','cm128_3'};
for i = 1:5
    fn = session_names{i};
    IDXtmp = IDX_final.(fn);
    IDXtmp(IDXtmp>46) = IDXtmp(IDXtmp>46) - 92;
    IDXtmp = IDXtmp + 46;
    IDX_final.(fn) = IDXtmp;
    
    
end

%% Eroding
fn = {'cm124_2','cm125_2','cm126_2','cm127_2','cm128_3'};
for i = 1:5
    IDXtmp = IDX_final.(fn{i});
    IDXnew = zeros(256,256);
    for j = 1:92
        ctmp = IDXtmp==j;
        ctmp = imfill(ctmp,'holes');
        ctmp = imerode(ctmp,ones(3));
        IDXnew(ctmp==1) = j;
    end
    IDX_final_new.(fn{i}) = IDXnew;
    
    
    
end
