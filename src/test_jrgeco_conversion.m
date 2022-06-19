function [Dr_best,Dg_best] = test_jrgeco_conversion(data,bl,Dr,Dg)

results = [];
for i = 1:numel(Dr)
    for j = 1:numel(Dg)
        tmp = (data.lime)./((data.red).^Dr(i)).*((data.green).^Dg(j));
        jrgeco{i,j} = tmp./repmat(nanmean(tmp(:,:,bl),3),[1 1 size(tmp,3)])-1;
        results(i,j) = nanstd(jrgeco{i,j}(:));
    end
end

[row,col] = find(results == min(min(results)));
Dr_best = Dr(row);
Dg_best = Dg(col);