% Compare centroids
s1 = {'State 1','State 2','State 3','State 4','State 5'};

s2 = {'onset','locomotion', 'offset', 'initial rest', 'sustained rest'};
load(fullfile('/local_mount/space/dingus/1/RS_analysis/Draft/4-dFC/results','resultsj'))
load(fullfile('/local_mount/space/dingus/1/RS_analysis/Draft/3-bFC/behavioral_states_new/states/10-s','states_jRGECO.mat'))
c1 = DFCopts.cntj;
c2 = reshape(finalStates(5).states,[92*92,5])';
c2 = c2([5 1 2 3 4],:);
figure;
suptitle('Comparison of state centroids: Behavior vs. dFC')
for i = 1:5
    subplot(2,5,i)
    imagesc(reshape(c1(i,:),[92 92]))
    axis image; axis off; caxis([0 1])
    title(s1{i})
    
    subplot(2,5,i+5)
    imagesc(reshape(c2(i,:),[92 92]))
    axis image; axis off; caxis([0 1])
    title(s2{i})
end
colormap jet