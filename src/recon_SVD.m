function dataout = recon_SVD(data,m)
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
rotf = getrotf(m);
m.baseline = baselinefromrot(rotf,100,1);
for d = 1:m.nLEDs
    dataout.(m.LEDs{d}) = reshape(data.C.(m.LEDs{d})*data.S.(m.LEDs{d}),...
        [m.height/m.dsf,m.width/m.dsf,size(data.S.(m.LEDs{d}),2)]);
    dataout.(m.LEDs{d}) = dataout.(m.LEDs{d}) - min(dataout.(m.LEDs{d})(:));
end
dataout.jrgeco = dataout.lime./((dataout.red.^m.Dr).*(dataout.green.^m.Dg));
[dataout.chbo,dataout.chbr,~] = convert_mariel_MIAO(dataout.green,dataout.red,'g','r',m.baseline,m.greenfilter);
bgGG = mean(dataout.jrgeco(:,:,m.baseline),3);
dataout.jrgeco = dataout.jrgeco./repmat(bgGG,[1 1 size(dataout.jrgeco,3)])-1;
dataout = structfun(@(x) (changem(x,[0 0],[NaN inf])),dataout,'UniformOutput',false);

