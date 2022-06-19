%% New IDX
addpath('/local_mount/space/chronos/2/RS_SVD_2020')
runnames = {'cm125_2_runD.mat','cm126_2_runD.mat','cm127_2_runD.mat','cm128_2_runD.mat'};
for i = 1:numel(runnames)
    clearvars -EXCEPT all_data all_m runnames i
    load(runnames{i})
    m.rotf = getrotf(m);
    m.baseline = baselinefromrot(m.rotf,100,1);
    for d = 1:m.nLEDs
        data.(m.LEDs{d}) = convertSVD(C.(m.LEDs{d}),S.(m.LEDs{d}),1:500,m.nanidx,256);
    end
    data.jrgeco = data.lime./((data.red.^m.Dr).*(data.green.^m.Dg));
    bgGG = mean(data.jrgeco(:,:,m.baseline),3);
    data.jrgeco = data.jrgeco./repmat(bgGG,[1 1 size(data.jrgeco,3)])-1;
    data.jrgeco = reshape(data.jrgeco,[256^2,size(data.jrgeco,3)]);
    jrgeco = data.jrgeco(:,1:11000);
    jrgeco(isinf(jrgeco)) = 0; jrgeco(isnan(jrgeco)) = 0;
    m.rotf = m.rotf(1:11000);
    clear data
    
    % hp filter
    jrgeco_filt = zeros(size(jrgeco))';
    fs = m.framerate/3;
    f_n = [0.5, 0.6]*2/fs; % stopband in normalized frame rate, so here it is 0.5Hz hp-filter
    a = [1 0];
    b_hp = firpm(500,[0 f_n 1], [0 0 1 1]); % filter order is 500
    non_nan_ind = ~isnan(jrgeco(:,1));
    jrgeco_filt(:,non_nan_ind) = filtfilt(b_hp,a,jrgeco(non_nan_ind,:)');
    jrgeco_filt = jrgeco_filt';
    jrgeco = jrgeco_filt;
    all_data.(m.mouse) = jrgeco;
    all_m.(m.mouse) = m;
end

%% Refine IDX
n = 4;
clear GUI_m
fn = fieldnames(all_m);
jrgeco = all_data.(fn{n});
m = all_m.(fn{n});


