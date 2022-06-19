function varargout = TF_GUI(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TF_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @TF_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function TF_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.IDX = varargin{1};
handles.refim = varargin{2};
axes(handles.axes1)
imagesc(handles.refim); axis image; axis off; colormap gray
im = getframe(gca);
handles.im = imresize(im.cdata,size(handles.IDX));
addlistener([handles.X,handles.Y,handles.R],'ContinuousValueChange',@(hObject, event) update_view(hObject,handles));
guidata(hObject, handles);
update_view(hObject,handles)
uiwait

function varargout = TF_GUI_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.TF;
close

function R_Callback(hObject, eventdata, handles)
update_view(hObject,handles)

function S_Callback(hObject, eventdata, handles)
update_view(hObject,handles)

function X_Callback(hObject, eventdata, handles)
update_view(hObject,handles)

function Y_Callback(hObject, eventdata, handles)
update_view(hObject,handles)

function update_view(hObject,handles)
TF = affine2d;
R = deg2rad(handles.R.Value);
X = handles.X.Value;
Y = handles.Y.Value;
S = str2num(handles.S.String);

% Rotate about center
RX = (128-128*cos(R)+128*sin(R));
RY = 128-128*sin(R)-128*cos(R);
TF.T = [S*cos(R)  sin(R)     0;
        -sin(R)   S*cos(R)   0;
        RX+X      RY-Y       1];

handles.TF = TF;
axes(handles.axes1); cla
handles.IDX_TF = imwarp(imerode(logical(imgradient_nic(handles.IDX)),ones(2)),handles.TF,'OutputView',imref2d(size(handles.IDX)));
imagesc(imoverlay(handles.im,handles.IDX_TF,'w')); axis image; axis off
guidata(hObject, handles);

function X_CreateFcn(hObject, eventdata, handles)

function Y_CreateFcn(hObject, eventdata, handles)

function S_CreateFcn(hObject, eventdata, handles)

function R_CreateFcn(hObject, eventdata, handles)

function closebutton_Callback(hObject, eventdata, handles)
update_view(hObject,handles)
uiresume
