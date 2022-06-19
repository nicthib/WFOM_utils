function stim_data = test_stim(m,data,Dr,Dg)
fn = fieldnames(data);

for i = 1:numel(m.stimframes)
    for j = 1:numel(fn)
        if i == 1
            stim_data.(fn{j}) = [];
            stim_data.(fn{j}) = data.(fn{j})(:,:,[-99:200]+m.stimframes(i));
        else
            stim_data.(fn{j}) = stim_data.(fn{j}) + data.(fn{j})(:,:,[-99:200]+m.stimframes(i));
        end
    end
end

% select stim region
imagesc(mean(stim_data.lime(:,:,100:130),3)-mean(stim_data.lime(:,:,1:70),3))
axis image; axis off; colormap gray; title('Pick response region')
whisk_roi = roipoly;

for i = 1:numel(Dr)
    for j = 1:numel(Dg)
        jrgeco = stim_data.lime./((stim_data.red.^Dr(i)).*(stim_data.green.^Dg(j)));
        jrgeco = jrgeco./repmat(mean(jrgeco(:,:,1:70),3),[1 1 size(jrgeco,3)])-1;
        whisk_TC = squeeze(nanmean(nanmean(jrgeco.*repmat(whisk_roi,[1 1 size(jrgeco,3)]))));
    end
end