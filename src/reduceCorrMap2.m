function corrMap_re = reduceCorrMap2(corrMap)
% reduceCorrMap reorder a FC map from an ROI-ROI map to a network-ROI
% input: corrMap (a 92x92 correlation matrix)
% output: corrMap_re (a 12x92 correlation matrix)

% Sub-groups defined by Nic
lM1 = 1:8;
rM1 = lM1 + 46;

lM2 = 9:16;
rM2 = lM2 + 46;

lS1 = 17:23;
rS1 = lS1 + 46;

lS2 = 24:27;
rS2 = lS2 + 46;

lS3 = 28:34;
rS3 = lS3 + 46;

lV = 35:46;
rV = lV + 46;


sets = {lM1, lM2, lS1, lS2, lS3, lV, rM1, rM2, rS1, rS2, rS3, rV};

corrMap_re = zeros(12,92);

% inter-network correlations
for i = 1:length(sets)
    for j = 1:92
        corrMap_re(i,j) = mean(corrMap(sets{i},j));
    end
end

end