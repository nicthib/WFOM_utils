figure
fn = fieldnames(results);
cmap = hsv(8); cmap(1,:) = 0;
c=1;
for i = 1:4
    cIDX{i} = imerode(logical(imgradient(results.(fn{i}).IDX_LR_refined)),ones(2));
    cIDX{i} = double(cIDX{i});
    cIDX_comps{i} = [];
    subplot(2,2,i)
    imagesc(cIDX{i})
    hold on
    axis image; axis off; colormap(cmap)
    
end

while true
    c = c+1;
    for i = 1:4
        subplot(2,2,i)
        [x,y] = getpts;
        x=round(x); y=round(y);
        cIDX_comps{i}{c-1} = [];
        for j = 1:numel(x)
            comp = results.(fn{i}).IDX_LR_refined(y(j),x(j));
            cIDX{i}(results.(fn{i}).IDX_LR_refined==comp) = c;
            cIDX_comps{i}{c-1}(end+1) = comp;
        end
        cla; imagesc(cIDX{i})
        title(numel(unique(cIDX_comps{i}{end})))
        caxis([0 8])
    end
end
%%
for i = 1:4
   cIDX{i} = cIDX{i}-1; 
end

%%
fn = fieldnames(results);
for i = 1:numel(fn)
    results.(fn{i}).IDX_final = zeros(256,256);
end
for i = 1:numel(cIDX_comps)
    for j = 1:numel(cIDX_comps{i})
        for k = 1:numel(cIDX_comps{i}{j})
            tmp = cIDX_comps{i}{j}(k);
            n = max(results.(fn{i}).IDX_final(:))+1;
            results.(fn{i}).IDX_final(results.(fn{i}).IDX_LR_refined==tmp) = n;
        end
    end
end

% Split LR
load('/local_mount/space/dingus/1/RS_analysis/registration/TF.mat','TF_w')
for i = 1:numel(fn)
    inc_idx = and(results.(fn{i}).IDX_final,results.(fn{i}).L)==1;
    results.(fn{i}).IDX_final(inc_idx==1) = results.(fn{i}).IDX_final(inc_idx==1) + max(results.(fn{i}).IDX_final(:));
    results.(fn{i}).IDX_final(results.(fn{i}).IDX_final==0) = NaN;
end

%%
mouse = 'cm128';
day = '8';
clear TF1 TF2
TF1 = TF_w.([mouse '_2_1']);
TF2 = TF_w.([mouse '_1_' day]);
TF1.tform.T = TF1.tform.T*TF2.tform.T
regd = imwarp_fr(results.([mouse '_2']).IDX_final,TF1);

close all
imagesc(logical(imgradient_nic(regd)))
axis image; axis off




