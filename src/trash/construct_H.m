function [H,Hstd] = construct_H(varargin)
IDX = varargin{1};
filename = varargin{2};
if numel(varargin) == 3
    savepath = varargin{3};
end
try
    H = []; Hstd = [];
    savefilename = [filename '_H.mat'];
    disp('Loading SVD components...')
    load(filename,'C','S','m');
    m.rotf = getrotf(m);
    m.baseline = baselinefromrot(m.rotf,100,1);
    disp(['Converting ' filename])
    disp('Converting SVD to raw LEDs...')
    for d = 1:m.nLEDs
        data.(m.LEDs{d}) = convertSVD(C.(m.LEDs{d}),S.(m.LEDs{d}),1:500,m.nanidx,256);
    end
    
    disp('Converting data...')
    data.jrgeco = data.lime./((data.red.^m.Dr).*(data.green.^m.Dg));
    [data.chbo,data.chbr,~] = convert_mariel_MIAO(data.green,data.red,'g','r',m.baseline,m.greenfilter);
    bgGG = mean(data.jrgeco(:,:,m.baseline),3);
    data.jrgeco = data.jrgeco./repmat(bgGG,[1 1 size(data.jrgeco,3)])-1;
    % reshape to speed up H
    data.jrgeco = reshape(data.jrgeco,[256^2,size(data.jrgeco,3)]);
    data.chbo = reshape(data.chbo,[256^2,size(data.chbo,3)]);
    data.chbr = reshape(data.chbr,[256^2,size(data.chbr,3)]);
    % GetHfromKmeans
    disp('Constructing H...')
    for e = 1:max(IDX(:))
        tmp = double(IDX == e);
        tmp(tmp == 0) = NaN; idx = ~isnan(reshape(tmp,[1 256^2]));
        H.jrgeco(e,:) = squeeze(nanmean(data.jrgeco(idx,:),1));
        H.chbo(e,:) = squeeze(nanmean(data.chbo(idx,:),1));
        H.chbr(e,:) = squeeze(nanmean(data.chbr(idx,:),1));
        Hstd.jrgeco(e,:) = squeeze(nanstd(data.jrgeco(idx,:),1));
        Hstd.chbo(e,:) = squeeze(nanstd(data.chbo(idx,:),[],1));
        Hstd.chbr(e,:) = squeeze(nanstd(data.chbr(idx,:),[],1));
    end
    m.DAQ = []; m.aux = []; % Free up space
    if numel(varargin) == 3
        save(fullfile(savedir,savefilename),'H','Hstd','m','-v7.3')
        disp(['Saved file ' savefilename ' in ' savedir])
    end
catch me
    disp(['Failed for ' filename, ', error is ' me.message])
end
        