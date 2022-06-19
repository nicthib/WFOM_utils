function regionalMatrix = group_corrs(FC,IDX,regions)
n = size(FC,1);
if FC(1,1) == 0
    FC = FC + eye(size(FC));
end
sz = size(IDX);
MM = zeros([sz,n]);
for i = 1:n
    tmp = zeros(sz);
    % Isolates each ROI from IDX map
    tmp(find(IDX==i))=1;
    % Put each ROI on own layer in 3D matrix
    MM(:,:,i) = tmp;
end
% regions = {[1:7], [8:13], [14:23], [24:27], [28:34], [35:46]... % left hemi
%     [47:53], [54:59], [60:69], [70:73], [74:80], [81:92]}; % right hemi
regionalMatrix = zeros(numel(regions),n);
for i = 1:length(regions)
    regionalMatrix(i,:) = mean(FC(regions{i},:)); % Mean across columns
end
