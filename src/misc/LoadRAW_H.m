function H = LoadData_RS_RAW(varargin)
m.fulldir = varargin{1};
[m.mouse,m.run,m.stim,m.CCDdir] = WFOM_splitdir(m.fulldir);
f = fieldnames(varargin{2});
for j = 1:length(f)
    m.(f{j}) = varargin{2}.(f{j});
end
m = ReadInfoFile(m);
try
    m = ReadAuxillary(m);
    disp('aux and DAQ data loaded')
catch
    disp('Unable to load aux and DAQ data')
end
m.rotf = getrotf(m);

%%%%%%%%%%%%%%%%%%%%%%%%%%% ZYLA LOAD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zylaInfoFilePath = fullfile(m.fulldir,'acquisitionmetadata.ini');
FID = fopen(zylaInfoFilePath, 'r');
zylaMetaData = fread(FID, '*char')';
fclose(FID);
AOIHeight_start = strfind(zylaMetaData, 'AOIHeight = ');
AOIWidth_start = strfind(zylaMetaData, 'AOIWidth = ');
AOIStride_start = strfind(zylaMetaData, 'AOIStride = ');
PixelEncoding_start = strfind(zylaMetaData, 'PixelEncoding = ');
ImagesPerFile_start = strfind(zylaMetaData, 'ImagesPerFile = ');
numDepths = str2double(zylaMetaData(AOIHeight_start+length('AOIHeight = '):...
    AOIWidth_start-1));
strideWidth = str2double(zylaMetaData(AOIStride_start+length('AOIStride = '):...
    PixelEncoding_start-1));
numLatPix = strideWidth/2;
numFramesPerSpool = str2double(zylaMetaData(ImagesPerFile_start+ length('ImagesPerFile = ')...
    :end));
numColumns = numDepths + m.offsetfactor;
numRows = numLatPix;    
for i = 1:length(dir(fullfile(m.fulldir,'*.dat')))
    temp = i-1; % added -1 here. this avoids not loading the 0000000000spool.dat file
    for j = 1:10
        a(i,j) = mod(temp, 10^j)/(10^(j-1));
        temp = temp-mod(temp, 10^j);
    end
    tempName = mat2str(a(i, :));
    tempName = tempName(2:end-1);
    tempName = tempName(find(tempName ~= ' '));
    tempName = [tempName 'spool.dat'];
    namesOut{i} = tempName;
end
if m.loadpct(1) == 0
    m.loadpct(1) = 1/numel(namesOut);
end
m.spoolsLoaded = round(numel(namesOut)*m.loadpct(1):round(numel(namesOut)*m.loadpct(2)));
filesToLoad = namesOut(m.spoolsLoaded);
FID = fopen(fullfile(regexprep(m.fulldir,'[_][0-9]$',''),'0000000000spool.dat'));
rawData = fread(FID, 'uint16=>uint16');
fclose(FID);
numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
try
    rawData = (reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
catch
    disp([sprintf('reshaping failed. Please check offset factor and adjust to proper value.\n Height = %i, Width = %i, factors = ',numRows,numColumns) mat2str(divisor(numel(rawData)/numRows))])
    return
end
% This switches height and width if nrot is odd
if mod(m.nrot,2) == 1
    height = m.height;
    m.height = m.width;
    m.width = height;
end
for i = 1:m.nLEDs
    F = rawData(:,:,i);
    F = F(1:m.height,1:m.width,:);
    if strcmp(m.bkgsub,'firstspool')
        m.bkg.([m.LEDs{i}]) = F;
    elseif ~isempty(m.bkgsub)
        m.bkg.([m.LEDs{i}]) = uint16(m.bkgsub);
    else
        m.bkg.([m.LEDs{i}]) = F*0;
    end
end
textprogressbar(sprintf(['Loading ' m.mouse ' ' m.run ' stim ' mat2str(m.stim) '\n']))
for i = 1:m.nLEDs  % preallocate
    data.(m.LEDs{i}) = (zeros(m.height/m.dsf,m.width/m.dsf,ceil(numFramesPerSpool.*length(filesToLoad)/m.nLEDs)));
end

for k = 1:m.nLEDs
    start = k;
    indexx.(m.LEDs{k}){1} = start:m.nLEDs:numFramesPerSpool;
    dsindex.(m.LEDs{k}){1} = 1:size(indexx.(m.LEDs{k}){1},2);
    indnum.(m.LEDs{k}){1} = length(dsindex.(m.LEDs{k}));
end

for i = 2:length(filesToLoad)
    for k = 1:m.nLEDs
        last = max(indexx.(m.LEDs{k}){i-1});
        skip = numFramesPerSpool-last;
        indexx.(m.LEDs{k}){i} = (m.nLEDs-skip-1)+[1:m.nLEDs:(numFramesPerSpool-(m.nLEDs-skip-1))]; % location in each batch of LEDs
        dsindex.(m.LEDs{k}){i} = max(dsindex.(m.LEDs{k}){i-1})+[1:size(indexx.(m.LEDs{k}){i},2)];
        indnum.(m.LEDs{k}){i} = length(dsindex.(m.LEDs{k}){i});
    end
end

for j = 1:length(filesToLoad)
    FID = fopen(fullfile(m.fulldir,filesToLoad{j}));
    rawData = fread(FID, 'uint16=>uint16');
    fclose(FID);
    numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
    numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
    if j == 1
        try
            rawData=(reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
        catch
            numColumns = numDepths + 1;
            numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
            numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
            rawData=(reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
        end
    else
        rawData=(reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
    end
    % Crop rawdata to original resolution
    rawData = rawData(1:m.height,1:m.width,:);
    for k = 1:m.nLEDs
        if m.dsf > 1
            data.(m.LEDs{k})(:,:,dsindex.(m.LEDs{k}){j}) = uint16(squeeze(mean(mean(reshape(rawData(1:m.height,1:m.width,indexx.(m.LEDs{k}){j}),[m.dsf,m.height/m.dsf,m.dsf,m.height/m.dsf,numel(indexx.(m.LEDs{k}){j})]),3),1))) - repmat(imresize(m.bkg.(m.LEDs{k}),1/m.dsf),[1 1 numel(indexx.(m.LEDs{k}){j})]);
        else
            data.(m.LEDs{k})(:,:,dsindex.(m.LEDs{k}){j}) = rawData(1:m.height,1:m.width,indexx.(m.LEDs{k}){j})-repmat(m.bkg.(m.LEDs{k}),[1 1 numel(indexx.(m.LEDs{k}){j})]);
        end
    end

    if mod(j,10)
        try
            textprogressbar(round(j*100/length(filesToLoad)));
        catch
            textprogressbar(sprintf(['Loading ' m.mouse ' ' m.run ' stim ' mat2str(m.stim) '\r']))
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotate
disp('Rotating LED data...')
for i = 1:m.nLEDs
    data.(m.LEDs{i}) = rot90(data.(m.LEDs{i}),m.nrot);
end
%

% LED reorder
if isfield(m,'LED_offset') && m.LED_offset ~= 0
    data_reorder = [];
    for i = 1:numel(m.LEDs)
       data_reorder.(m.LEDs{i}) = data.(m.LEDs{mod(i-1+m.LED_offset,m.nLEDs)+1});
    end
    data = data_reorder;
    clear data_reorder
end

% IDX registration
m.FIXED = data.green(:,:,1);
m.MOVING = m.regim;
m.TF = register_WFOM_session(m.FIXED,m.MOVING);
m.IDX = roundIDX(imwarp(m.IDX,m.TF,'OutputView',imref2d(size(m.MOVING))));
%

% H extraction and spectra check
disp('Extracting H...')
fn = fieldnames(data);
for i = 1:numel(fn)
    H.(fn{i}) = getHfromKmeans(data.(fn{i}),m.IDX,0);
end
%