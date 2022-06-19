function [data,m] = jrgeco_correction(data,m)
arguments
    data struct
    m struct
end
F = data.lime./((data.red.^m.Dr).*(data.green.^m.Dg));
fn = fieldnames(data);

if size(F) == 3
    for i = 1:numel(fn)
        m.bl_im.(fn{i}) = mean(data.(fn{i})(:,:,m.baseline),3);
    end
    bl = nanmean(F(:,:,m.baseline),3);
    dF = F./repmat(bl,[1 1 size(F,3)]); % F/F0
else
    bl = nanmean(F(:,m.baseline),2);
    dF = F./repmat(bl,[1 size(F,2)]); % F/F0
end

data.jrgeco = dF-1;
F = F/max(F(:)); % raw F normalized
F_grad = zeros(size(F));
dF_grad = zeros(size(F));
if size(F)==3
    for k = 1:size(dF,3)
        dF_grad(:,:,k) = imgradient(dF(:,:,k));
        F_grad(:,:,k) = imgradient(F(:,:,k));
    end
    sz = size(data.(fn{1}),1); % True resolution
    sz_20 = round(sz/5); % 20% of resolution
    grad_ran = sz_20:sz-sz_20;
    dF_grad = dF_grad(grad_ran,grad_ran,:); % F/F0 gradient
    F_grad = F_grad(grad_ran,grad_ran,:); % Raw F gradient
    m.jrgeco_quality_TC = squeeze(nanmean(nanmean(F_grad)));
    m.jrgeco_quality_im = max(cat(3,dF_grad,F_grad),[],3);
    m.jrgeco_conversion_score = mean(m.jrgeco_quality_im(:));
end


