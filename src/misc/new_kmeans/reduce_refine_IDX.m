function IDX_LR_refined = reduce_refine_IDX(IDX_L,IDX_R,IDX_LR)

figure('Position',[0 0 1000 500])
IDX_LR_refined = zeros(256,256);
currcomp = zeros(256,256);
i = min(find(~ismember(1:50,unique(IDX_LR_refined(:)))));
subplot(121); imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
title('Press Up arrow to store component.')
subplot(122); imagesc(IDX_LR_refined); axis image; axis off
while true
    try
        [x,y] = ginput(1); x = round(x); y = round(y);
        if IDX_L(y,x) >= 1
            currcomp(IDX_L==IDX_L(y,x)) = ~currcomp(y,x)*2;
        else
            currcomp(IDX_R==IDX_R(y,x)) = ~currcomp(y,x)*2;
        end
        subplot(121); imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
        
%         subplot(132);
%         if sum(currcomp(:)) > 0
%             c = corr(H',getHfromKmeans(data,currcomp/2,0)');
%             cIDX = zeros(256,256);
%             for j = 1:max(IDX_LR(:))
%                 cIDX(IDX_LR==j) = c(j);
%             end
%             imagesc(cIDX); axis image; axis off; caxis([0 1])
%         end
        title('Press Up arrow to store component.')
        k = waitforbuttonpress;
        v = double(get(gcf,'CurrentCharacter'));
        if v == 30 % U arrow (STORE)
            IDX_LR_refined = IDX_LR_refined + (currcomp/2).*i;
            subplot(122); imagesc(IDX_LR_refined); axis image; axis off
            currcomp = zeros(256,256);
            i = min(find(~ismember(1:50,unique(IDX_LR_refined(:)))));
        elseif v == 31 % D arrow (DELETE)
            IDX_LR_refined(find(currcomp)) = 0;
            currcomp = zeros(256,256);
            subplot(122); imagesc(IDX_LR_refined); axis image; axis off
            % Update
            subplot(121); imshowpair(currcomp+logical(IDX_LR_refined),logical(imgradient(IDX_LR))); axis image; axis off
            i = min(find(~ismember(1:50,unique(IDX_LR_refined(:)))));
        end
    catch
        break
    end
end