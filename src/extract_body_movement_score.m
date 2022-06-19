function [B_TC] = extract_body_movement_score(runname)
%runname = 'cm125_4_runC';
filepath = '/local_mount/space/dingus/1/RS_analysis/analysis_dec2020/H';
DLCdir = '/local_mount/space/dingus/1/RS_analysis/DLC/final';
load(fullfile(filepath,[runname '_H.mat']),'m')
B = extract_DLC_vars(DLCdir,runname);
B_H = [];
B_H.x = [];
B_H.y = [];
B_H.label = {};
fn = fieldnames(B);
for i = 1:numel(fn)
    body_part = B.(fn{i});
    if ~isempty(body_part) && isstruct(body_part)
        fn2 = fieldnames(body_part);
        for j = 1:numel(fn2)
            tmpvar = body_part.(fn2{j});
            if size(tmpvar,2) == 3
                B_H.x(end+1,:) = tmpvar(:,1)';
                B_H.y(end+1,:) = tmpvar(:,2)';
                B_H.x(end,tmpvar(:,3)<.95) = NaN;
                B_H.y(end,tmpvar(:,3)<.95) = NaN;
                B_H.label{end+1,1} = fn2{j};
            end
        end
    end
end

B_full = nanzscore([B_H.x; B_H.y]')';
B_std = movstd(B_full',20)';
B_TC = nanmean(B_std(:,1:3:end));
%close all
%figure
%plot(B_TC)
%hold on
%plot(abs(m.rotf/max(m.rotf)))
%legend('Body movement score','Locomotion')
