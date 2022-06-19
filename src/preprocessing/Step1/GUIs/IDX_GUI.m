function varargout = IDX_GUI(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IDX_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @IDX_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function IDX_GUI_OpeningFcn(hObject, ~, h, varargin)
h.output = hObject;
warning('off','stats:kmeans:MissingDataRemoved')

% Inputs
h.m = varargin{1};
h.jrgeco = varargin{2};
h.m.sz = h.m.height/h.m.dsf;

% Defaults
if numel(size(h.jrgeco)) == 3
    h.jrgeco = reshape(h.jrgeco,[size(h.jrgeco,1)*size(h.jrgeco,2) size(h.jrgeco,3)]);
end 
h.IDX = ones(h.m.sz,h.m.sz);
h.BW = ones(h.m.sz,h.m.sz);
h.BW_L = logical([ones(h.m.sz,h.m.sz/2) zeros(h.m.sz,h.m.sz/2)]);
h.BW_R = ~h.BW_L;
h.L = h.BW_L;
h.R = h.BW_R;
h.kmeans_idx = [0 0];
h.midptx = [h.m.sz/2 h.m.sz/2];
h.midpty = [h.m.sz 0];
h.IDX_LR_refined = h.IDX;
h.fr = h.m.framerate/h.m.nLEDs;
h.t = linspace(0,numel(h.m.rotf)/h.fr,numel(h.m.rotf));
h.TF = affine2d;
h.refim = h.m.refim;

guidata(hObject, h);
update_plots(h);
uiwait

function varargout = IDX_GUI_OutputFcn(hObject, ~, h)
if size(h.IDX,3) == 1
    IDX = sortkmeans(reduceIDX(h.IDX_LR_refined)); % Remove gaps in ROI values and sort top to bottom 
    h.IDX_LR_refined = IDX + max(IDX(:))*and(IDX,~h.L); % Bilateralize
end
fn = {'IDX_LR_refined','BW','BW_L','BW_R','L','R'};
m = h.m;
for i = 1:numel(fn)
    m.(fn{i}) = roundIDX(imwarp(h.(fn{i}),h.TF,'OutputView',...
    imref2d([h.m.sz,h.m.sz]))); % Apply straightening
end
m.refim = imwarp(h.refim,h.TF,'OutputView',...
    imref2d([h.m.sz,h.m.sz])); % Apply straightening
m.IDX = zeros(size(h.IDX));
for i = 1:size(h.IDX,3)
    m.IDX(:,:,i) = roundIDX(imwarp(h.IDX(:,:,i),h.TF,'OutputView',...
    imref2d([h.m.sz,h.m.sz]))); % Apply straightening
end
varargout{1} = m;
closereq

function run_cluster_Callback(hObject, ~, h)
h.BW_L = and(and(h.BW,h.IDX(:,:,1)),h.L);
h.BW_R = and(and(h.BW,h.IDX(:,:,1)),h.R);
h.exemplar1 = h.jrgeco(:,h.kmeans_idx(1):h.kmeans_idx(2));
if h.show_L.Value
    h.exemplar1(h.BW_L==0,:) = NaN;
elseif h.show_R.Value
    h.exemplar1(h.BW_R==0,:) = NaN;
else
    h.exemplar1(h.BW==0) = NaN;
end

k_num = str2num(h.n_clust.String);

for i = 1:numel(k_num)
    if numel(k_num) == 1
        missing_K = str2num(h.n_clust.String) - numel(unique(h.IDX_LR_refined(:))) + 1;
        if h.show_LR.Value
            missing_K = missing_K*2;
        end
    else
        missing_K = k_num(i) + k_num(i)*h.show_LR.Value;
    end
    IDXtmp = kmeans(h.exemplar1,missing_K,'Distance','correlation');
    IDXtmp(isnan(IDXtmp)) = 0;
    h.IDX(:,:,i) = reshape(IDXtmp,[h.m.sz, h.m.sz]);

    if h.show_L.Value || h.show_R.Value
        % Get exemplar again
        H = getHfromKmeans(h.exemplar1,h.IDX(:,:,i),1);
        h.exemplar2 = h.jrgeco(:,h.kmeans_idx(1):h.kmeans_idx(2));
        if h.show_L.Value
            h.exemplar2(h.BW_R==0,:) = NaN;
        elseif h.show_R.Value
            h.exemplar2(h.BW_L==0,:) = NaN;
        end
        
        IDX_i = kmeans(h.exemplar2,missing_K,'Distance','correlation','Start',H);
        IDX_i = reshape(IDX_i,[h.m.sz h.m.sz]);
        IDX_i(isnan(IDX_i)) = 0;
        h.IDX(:,:,i) = h.IDX(:,:,i)+IDX_i;
    end
end
update_plots(h)
guidata(hObject, h);

function select_kmeans_idx_Callback(hObject, ~, h)
axes(h.rot_ax)
[h.kmeans_idx,~] = ginput(2);
h.kmeans_idx = round(h.kmeans_idx*h.fr);
update_plots(h)
guidata(hObject, h);

function pick_midline_Callback(hObject, ~, h)
axes(h.IDX_ax)
[h.midptx,h.midpty] = getpts;
if numel(h.midptx) ~=2
    errordlg('Please pick only 2 points!')
    return
end

% re-evaluate line at top and bottom of frame
ri = diff(h.midpty);
ru = diff(h.midptx);
sl = ri/ru;
b = h.midpty(1)-sl*h.midptx(1);
h.midpty = [0 h.m.sz];
h.midptx = ([0 h.m.sz] - b)/sl;
h.L = poly2mask([h.midptx 0 0],[h.midpty h.m.sz 0],h.m.sz,h.m.sz);
h.R = ~h.L;
update_plots(h)
guidata(hObject, h);

function update_plots(h)
axes(h.rot_ax); cla;
plot(h.t,h.m.rotf); hold on; yticks([])
line([h.kmeans_idx(1) h.kmeans_idx(1)]/h.fr,[min(h.m.rotf) max(h.m.rotf)],'LineStyle','--','Color',[0 0 0])
line([h.kmeans_idx(2) h.kmeans_idx(2)]/h.fr,[min(h.m.rotf) max(h.m.rotf)],'LineStyle','--','Color',[0 0 0])
xlabel('Time (sec)')

axes(h.IDX_ax); cla
if h.show_L.Value
    tmpBW = h.BW.*h.L;
elseif h.show_R.Value
    tmpBW = h.BW.*h.R;
else
    tmpBW = h.BW.*or(h.R,h.L);
end
imagesc(h.IDX(:,:,1).*tmpBW); caxis([0 max(max(h.IDX(:,:,1)))]); hold on;
line(h.midptx,h.midpty,'LineStyle','--','Color',[1 1 1])
colormap jet; axis image; axis off
title('IDX')

IDXtmp = roundIDX(imwarp(h.IDX_LR_refined,h.TF,'OutputView',...
    imref2d(size(h.IDX_LR_refined))));
axes(h.IDX_out); imagesc(IDXtmp);
axis image; axis off; colormap jet

function show_L_Callback(hObject, ~, h)
update_plots(h)

function show_R_Callback(hObject, ~, h)
update_plots(h)

function show_LR_Callback(hObject, ~, h)
update_plots(h)

function crop_IDX_Callback(hObject, ~, h)
a=figure('Position',[200 200 700 700]);
imagesc(h.IDX(:,:,1)); axis image; colormap jet; axis image
roi = drawassisted;
h.BW = createMask(roi);
close(a)
guidata(hObject, h);
update_plots(h)

function quick_crop_Callback(hObject, ~, h)
a=figure;
imagesc(h.IDX(:,:,1)); axis image; colormap jet; axis image
roi = drawassisted;
BW = createMask(roi);
h.BW(BW==1) = 0;
close(a)
guidata(hObject, h);
update_plots(h)

function assign_clusts_Callback(hObject, ~, h)
if size(h.IDX,3) > 1
    edrrordlg('You can''t do this with more than 1 cluster number!')
    return
end
IDX_LR_refined = zeros(h.m.sz,h.m.sz);
IDX_LR = h.IDX; IDX_LR(isnan(IDX_LR)) = 0;
IDX_L = h.IDX; IDX_L(h.BW_L==0) = 0;
IDX_R = h.IDX; IDX_R(h.BW_R==0) = 0;
BW_L = h.BW_L;
BW_R = h.BW_R;
currcomp = zeros(h.m.sz,h.m.sz);
comp_candidates = 1:str2num(h.n_clust.String);
i = min(find(~ismember(comp_candidates,unique(IDX_LR_refined(:)))));

while true
    singlecomp = zeros(h.m.sz,h.m.sz);
    axes(h.IDX_ax); cla
    imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
    [xi,yi] = getpts; xi = round(xi); yi = round(yi);
    if isempty(xi)
        break
    end
    for j = 1:numel(xi)
        x = xi(j); y = yi(j);
        if BW_L(y,x) == 1
            singlecomp(IDX_L==IDX_L(y,x)) = 1;
            currcomp(IDX_L==IDX_L(y,x)) = ~currcomp(y,x)*2;
        elseif BW_R(y,x) == 1
            singlecomp(IDX_R==IDX_R(y,x)) = 1;
            currcomp(IDX_R==IDX_R(y,x)) = ~currcomp(y,x)*2;
        end
    end
    imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
    
    if sum(sum(and(currcomp,IDX_LR_refined))) == 0 
        IDX_LR_refined = IDX_LR_refined + (currcomp/2).*i;
        currcomp = zeros(h.m.sz,h.m.sz);
        axes(h.IDX_ax); imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
        axes(h.IDX_out); imagesc(IDX_LR_refined); axis image; axis off
        i = min(find(~ismember(comp_candidates,unique(IDX_LR_refined(:)))));
    else
        IDX_LR_refined(find(currcomp)) = 0;
        currcomp = zeros(h.m.sz,h.m.sz);
        axes(h.IDX_ax); imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
        axes(h.IDX_out); imagesc(IDX_LR_refined); axis image; axis off
        i = min(find(~ismember(comp_candidates,unique(IDX_LR_refined(:)))));
    end
end
if sum(h.IDX_LR_refined(:)) > 0
    IDX_LR_refined(IDX_LR_refined~=0) = IDX_LR_refined(IDX_LR_refined~=0)+max(h.IDX_LR_refined(:));
end
h.IDX_LR_refined = h.IDX_LR_refined+IDX_LR_refined;
h.IDX(IDX_LR_refined > 0) = 0;

guidata(hObject, h);
update_plots(h)

function straighten_Callback(hObject, ~, h)
bkg = zeros(size(h.IDX_LR_refined));
bkg(1:end,size(bkg,2)/2-1:size(bkg,2)/2+1) = 1;
if size(h.IDX,3) == 1
    h.TF = TF_GUI(h.IDX_LR_refined,bkg); 
else
    h.TF = TF_GUI(h.IDX(:,:,1),bkg); 
end
guidata(hObject, h);
update_plots(h)

function save_IDX_Callback(hObject, ~, h)
update_plots(h);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n_clust_Callback(hObject, ~, h)

function n_clust_CreateFcn(hObject, ~, h)

% function getSTD_Callback(hObject, ~, h)
% axes(h.IDX_out)
% silhL = IDX_H_std(h.IDX_LR_refined.*h.L,h.jrgeco(:,h.kmeans_idx(1):h.kmeans_idx(2)));
% silhR = IDX_H_std(h.IDX_LR_refined.*h.R,h.jrgeco(:,h.kmeans_idx(1):h.kmeans_idx(2)));
% silhR(isnan(silhR)) = 0; silhL(isnan(silhL)) = 0;
% h.silh = silhL+silhR;
% imagesc(h.silh)
% axis image; axis off; colorbar; caxis([0 1])
% 
% axes(h.IDX_ax)
% imagesc(h.IDX_LR_refined)
% caxis([min(h.IDX_LR_refined(:)) max(h.IDX_LR_refined(:))])
% guidata(hObject, h);
