function behavior = get_behavior(fn)
behavior = []; fn_DLC = strrep(fn,'_H.mat','');
try load('/local_mount/space/dingus/1/RS_analysis/DLC/final/pupil_032521.mat',fn_DLC)
    [ptmp,behavior.pupil_perf] = convertpupilcoords(eval(fn_DLC),.95);
     disp(['Pupil performance: ' mat2str(round(perf,2))])
    behavior.pupil = sepblockfun(ptmp(1:11980*3,3),[3 1],'mean');
catch; end

try load('/local_mount/space/dingus/1/RS_analysis/DLC/final/whisker1_060721.mat',fn_DLC)
    tmp = eval(fn_DLC);
    wn = {'most_prominent_whisker-1','most_prominent_whisker-2','most_prominent_whisker-3'};
    for f = 1:numel(wn)
        tmp2 = tmp.(wn{f});
        if ndims(tmp2) == 3
            tmp2 = squeeze(tmp2(1,:,:));
        end
        %dx = diff(tmp2(:,1));
        %dy = diff(tmp2(:,2));
        w_p = sepblockfun(tmp2(1:11980*3,3),[3 1],'mean');
        %w = sqrt(dx.^2+dy.^2);
        w = tmp2(:,1);
        whiskers(f,:) = sepblockfun(w(1:11980*3),[3 1],'mean');
    end
    behavior.whisk = nanmean(whiskers,1);
    behavior.whisk_p = w_p;
catch; end

%load(fn,'m')
%behavior.rotf = smooth(m.rotf',20);