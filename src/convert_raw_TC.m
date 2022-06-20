% Converts a TC using given metadata. Designed for RS project.
function dataout = convert_raw_TC(data,m,showsummary)
%%LP Filter
zlen = 1000; n = size(data.red,1);
data.red = [zeros(n,zlen),data.red,zeros(n,zlen)];
data.green = [zeros(n,zlen),data.green,zeros(n,zlen)];

m.lpFilt = designfilt('lowpassiir','FilterOrder',20,...
'PassbandFrequency',2,'PassbandRipple',0.2,'SampleRate',19.9668);

for n = 1:size(data.red,1)
    data.red(n,:) = filtfilt(m.lpFilt,data.red(n,:));
    data.green(n,:) = filtfilt(m.lpFilt,data.green(n,:));
end

data.red = data.red(:,zlen+1:end-zlen);
data.green = data.green(:,zlen+1:end-zlen);

% get bkg, then convert to bkg vals
bkgH = getHfromKmeans(double(repmat(imresize(m.bkg.lime,.5),[1 1 2])),m.IDX,0);
bkgH = bkgH(:,1);

% subtract bkg from raw TCs
fn = {'lime','green','red'};
for i = 1:numel(fn)
    data.(fn{i}) = data.(fn{i}) - repmat(bkgH,[1 size(data.(fn{i}),2)]);
end

% get dynamic baseline
dbl={}; idx = 1;
dbl.red = [];
dbl.green = [];
dbl.lime  = [];
dbl.idx = [];
blW = 50;
blB = 50;
blR = 500;
if numel(m.rotf) < 100
    m.rotf = getrotf(m);
end
for i=1:blR:numel(m.rotf)-blR
    bltestidx = i:i+blR;
    blidx = baselinefromrot(m.rotf(bltestidx),blW+2*blB,1)+i-1;
    if ~isempty(blidx)
        blidx = blidx(blB:end-blB-1);
        dbl.lime(:,idx) = mean(data.lime(:,blidx),2);
        dbl.green(:,idx) = mean(data.green(:,blidx),2);
        dbl.red(:,idx) = mean(data.red(:,blidx),2);
        dbl.idx(:,idx) = i+blB;
        idx = idx+1;
    end
end
% manual extrap with nearest neighbor
% first frame
% dbl.lime(:,end+1) = dbl.lime(:,1);
% dbl.green(:,end+1) = dbl.green(:,1);
% dbl.red(:,end+1) = dbl.red(:,1);

if numel(dbl.idx) < 2
    disp('Not enough baseline candidates found!')
    dataout = [];
    return
end

% apply dynamic baseline
dbl_full={};
dbl_full.red = [];
dbl_full.green = [];
dbl_full.lime  = [];
for i = 1:size(data.red,1)
    dbl_full.red(i,:) = interp1(dbl.idx,dbl.red(i,:),1:size(data.red,2),'linear','extrap');
    dbl_full.green(i,:) = interp1(dbl.idx,dbl.green(i,:),1:size(data.red,2),'linear','extrap');
    dbl_full.lime(i,:) = interp1(dbl.idx,dbl.lime(i,:),1:size(data.red,2),'linear','extrap');
end

% files for hb conversion
m.spectra_root = '/local_mount/space/juno/1/Software/MIAO/Spectra';
files = {};
files.spectra_file = m.spectra_file;
files.Hb_file = fullfile(m.spectra_root,'Hb_spectra.mat');
files.dpffs = fullfile(m.spectra_root,'dpffsMW_400to700nm.mat');

dataout = {};

% Hb Conversion
[dataout.chbo,dataout.chbr,~] = convert_Hb_2D(data.green,data.red,'g','r',dbl_full.green,dbl_full.red,files);
%

% jRGECO conversion
dataout.jrgeco = jrgeco_correction_2D(data,m,dbl_full);

% Save dynamic baseline
dataout.dbl = dbl;

if showsummary
    figure
    subplot(211)
    plot(data.green'-repmat(data.green(:,100),[1 size(data.green,2)])');

    subplot(212)
    plot(m.rotf)
end