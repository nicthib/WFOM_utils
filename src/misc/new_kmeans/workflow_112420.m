% Extract running exemplars

addpath('/local_mount/space/chronos/2/RS_SVD_2020')
mouse = 'cm125_2';
runname = [mouse '_runD.mat'];
% Get movement epochs
load(runname,'m')
rotf = getrotf(m);
rotf(rotf<0) = 0;
moveidx = find(smooth(rotf,50)>1);
epochs = {}; epoch = [];
while numel(moveidx) > 1
    d = diff(moveidx(1:2));
    if d == 1
        epoch = [epoch moveidx(1)];
    else
        epochs{end+1} = epoch;
        epoch = [];
    end
    moveidx(1) = [];
end
epochs(find(cellfun(@numel,epochs) < 100)) = [];
%epochs(find(cellfun(@numel,epochs) > 400)) = [];

AS_epochs.(mouse) = epochs;

%%
load(runname)
rotf = getrotf(m);
m.baseline = baselinefromrot(rotf,100,.5);
for d = 1:m.nLEDs
    data.(m.LEDs{d}) = convertSVD(C.(m.LEDs{d}),S.(m.LEDs{d}),1:500,m.nanidx,256);
end
data.jrgeco = data.lime./((data.red.^m.Dr).*(data.green.^m.Dg));
bgGG = mean(data.jrgeco(:,:,m.baseline),3);
data.jrgeco = data.jrgeco./repmat(bgGG,[1 1 size(data.jrgeco,3)])-1;
data.jrgeco = reshape(data.jrgeco,[256^2,size(data.jrgeco,3)]);
jrgeco = data.jrgeco(:,1:11000);
jrgeco(isinf(jrgeco)) = 0; jrgeco(isnan(jrgeco)) = 0;
rotf = rotf(1:11000);
clear data

%% Try clustering again
load('IDX_refined.mat')
load('IDX_final.mat','IDX_final')
fn = 'cm126_2';
BW_LR = bwareaopen(and(~logical(IDX_refined.(fn)),logical(IDX_final.(fn).LR)),1000);
%exemplar_cl = jrgeco(:,5416:6719);
exemplar_cl = exemplars.(fn).RS;
exemplar_cl(BW_LR(:)~=1,:) = NaN;
missing_K = 50-max(IDX_refined.(fn)(:));
IDX_somat = kmeans(exemplar_cl,missing_K*2,'Distance','correlation');
IDX_somat = reshape(IDX_somat,[256 256]);
close all; imagesc(IDX_somat)

%%
IDX_LR = IDX_somat; IDX_LR(isnan(IDX_LR)) = 0;
IDX_L = IDX_LR.*IDX_final.(fn).BW_L;
IDX_R = IDX_LR.*IDX_final.(fn).BW_R;
H = getHfromKmeans(exemplar_cl,IDX_LR,0);
reduce_refine_IDX(IDX_L,IDX_R,IDX_LR,H,exemplar_cl)

