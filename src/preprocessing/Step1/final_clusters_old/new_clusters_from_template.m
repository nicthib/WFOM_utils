function m = new_clusters_from_template(session_name,run_name,m,ref_data)

% mt (mouse target) is always cm125
% ms is the mouse of interest
IDX_start = ref_data.IDX_final;
p = [sin(linspace(0,2*pi,9))' cos(linspace(0,2*pi,9))']*120+120;
p=p(1:8,:);

load_dir = findmousefolder(session_name,run_name,1);
[~,ms] = LoadData_v21(load_dir,m);
tmp = max(ms.green,[],3); ms.refim = uint8(round(255*tmp/max(tmp(:))));
mt.p = p; ms.p = p;
mt.refim = uint8(round(255*ref_data.refim/max(ref_data.refim(:))));

% Select registration points
% Here, we show the "fiducial" points on cm125's mask, and you are expected
% to adjust the othe fiducials on the right image to match. We use this to
% perform a cross mouse registration transformation matrix.
[mt.p,ms.p] = cpselect(mt.refim,ms.refim,mt.p,ms.p,'Wait',true);

% Convert fiduciary points to masks
idx = convhull(mt.p(:,2),mt.p(:,1));
mt.regmask = uint8(poly2mask(mt.p(idx,1),mt.p(idx,2),256,256));
idx = convhull(ms.p(:,2),ms.p(:,1));
ms.regmask = uint8(poly2mask(ms.p(idx,1),ms.p(idx,2),256,256));

% Register masks
[optimizer, metric] = imregconfig('monomodal');
TF = imregister_nic(mt.regmask,ms.regmask,'affine',optimizer,metric);

% subplot(121)
% imshowpair(mt.refim,ms.refim)
% title('Before')
% 
% subplot(122)
% imshowpair(imwarp_fr(mt.refim,TF),ms.refim)
% title('After')
% a = questdlg('Is this good?','yes','no');
% if strcmp(a,'no')
%     return
% end

m.IDX_final = roundIDX(imwarp_fr(IDX_start,TF));
m.refim = mt.refim;
m.p = ms.p;

% m = m_default; m.outputs = 'gn'; %m.noload = 1;
% m.loadpct = [0 1]; m.BW = ones(256,256);
% [m,data] = LoadData_v2(load_dir,m);
% m.aux = []; m.DAQ = [];
% rotf = smooth(abs(m.rotf),200);
% m.score1 = []; m.score2 = [];  m.score3 = []; 
% m.epochs = [];
% 
% % remove epochs with running
% epoch_length = 1200;
% for i = 1:200:10000
%     if max(rotf(i:i+1200)) < .5
%        m.epochs(end+1) = i;
%     end
% end
% m.subregions = [9 16 9 5 7];
% m.refim = data.green(:,:,1);
% 
% for i = 1:numel(m.epochs)
%     m.kmeans_idx = m.epochs(i):m.epochs(i)+epoch_length;
%     jrgeco = data.jrgeco(:,:,m.kmeans_idx);
%     IDX_end(IDX_end==0) = NaN;
%     BW = double(IDX_end>0); BW(BW==0) = NaN;
%     jrgeco = jrgeco.*repmat(BW,[1 1 size(jrgeco,3)]);
%     m.H = getHfromKmeans(jrgeco,IDX_end,1);
%     m.IDX(:,:,i) = reshape(kmeans(reshape(jrgeco,[256^2 size(jrgeco,3)]),size(m.H,1),'Distance','correlation','Start',m.H),[256 256]);
%     IDX2 = m.IDX(:,:,i);
%     
%     % Metrics
%     % Check overlap
%     denom = sum(IDX2(:)>=1);
%     numer = sum(abs(IDX_end(:)-IDX2(:))==0);
%     m.score1(i) = round(numer/denom,2);
%     
%     % Check graph similarity
%     IDX1_n = IDXtodigraph(IDX_end);
%     IDX2_n = IDXtodigraph(IDX2);
%     m.score2(i) = round(graph_similarity(IDX1_n,IDX2_n),3);
% 
%     % Check centroid nearness
%     m.score3(i) = round(centroid_nearness(IDX_end,IDX2),1);
%     
%     % Display scores
%     disp(['Score1: ' mat2str(m.score1(i)) ', Score2: ' mat2str(m.score2(i)) ', Score3: ' mat2str(m.score3(i))])
% end
% IDX = m.IDX_best;
% m.IDX_best = m.IDX(:,:,find(m.score2==max(m.score2),1));
% 
% close all
% imagesc(logical(imgradient_nic(m.IDX_best))); axis image; axis off
% title('Draw a boundary around the left hemisphere')
% BW_L = roipoly;
% IDX(IDX>46) = IDX(IDX>46) - 46;
% m.IDX_final = IDX.*BW_L+((IDX+46).*~BW_L);
% 
% subplot(121)
% showIDX_labels_subregions(IDX_start,max(IDX_start(:)),[0 0 1]);
% title('Template')
% 
% subplot(122)
% showIDX_labels_subregions(m.IDX_final,max(m.IDX_final(:)),[0 0 1])
% title('Result')