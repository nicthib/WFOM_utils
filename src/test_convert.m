function Q = test_convert(data,D,m,metric)
if strcmp(metric,'RS')
    tmp = data.lime./((data.green.^D(2)).*(data.red.^D(1)));
    bl = nanmean(tmp(:,:,m.baseline),3);
    F_F0 = tmp./repmat(bl,[1 1 size(tmp,3)]); % F/F0
    F = tmp/max(tmp(:)); % raw F normalized
    F_F0 = F_F0(:,:,1:1000);
    F = F(:,:,1:1000);
    F_F0_grad = zeros(size(F_F0));
    F_grad = zeros(size(F));
    for k = 1:size(F_F0,3)
        F_F0_grad(:,:,k) = imgradient(F_F0(:,:,k));
        F_grad(:,:,k) = imgradient(F(:,:,k));
    end
    grad = max(cat(3,F_F0_grad,F_grad),[],3);
    Q = sum(grad(:)>.3)/numel(grad(:));
    
    % Show results
    figure(1)
    subplot(211)
    imagesc(grad); title(['Dr = ' mat2str(round(D(1),2)) ', Dg = ' mat2str(round(D(2),2)) ', Q = ' mat2str(round(Q,2))])
    axis image; axis off; colormap gray; colorbar; caxis([0 1])
    subplot(212)
    hold on;
    cmap = parula(100);
    a = gca;
    scatter(D(1),D(2),'filled','MarkerFaceColor',cmap(mod(numel(a.Children),100)+1,:))
    xlim([0 2])
    ylim([0 2])
    xlabel('Dr')
    ylabel('Dg')
    
elseif strcmp(metric,'STIM')
    jrgeco = m.stim_data.lime./((m.stim_data.red.^D(1)).*(m.stim_data.green.^D(2)));
    jrgeco = jrgeco./repmat(mean(jrgeco(:,:,1:99),3),[1 1 size(jrgeco,3)])-1;
    %whisk_TC = squeeze(nanmean(nanmean(jrgeco.*repmat(m.whisk_roi,[1 1 size(jrgeco,3)]))));
    whisk_TC = squeeze(nanstd(nanstd(abs(jrgeco))));
    Q = mean(abs(whisk_TC(130:200)))/max(whisk_TC);    
    %Q1 = mean(abs(whisk_TC(200:300)))/max(whisk_TC);
    %Q2 = -min(0,min(whisk_TC(100:300)))/max(whisk_TC);
    %Q = mean([Q1 Q2]);
    % baseline score is defined as abs(M_pre-M_post)/max(TC)
    % Best score is 0, theoretical worst is 1
    figure(1)
    
    subplot(211)
    plot(linspace(-5,10,300),whisk_TC)
    title(['Dr = ' mat2str(round(D(1),2)) ', Dg = ' mat2str(round(D(2),2)) ', Q = ' mat2str(round(Q,2))])
    
    subplot(212)
    hold on;
    cmap = parula(100);
    a = gca;
    scatter(D(1),D(2),'filled','MarkerFaceColor',cmap(mod(numel(a.Children),100)+1,:))
    xlim([0 2])
    ylim([0 2])
    xlabel('Dr')
    ylabel('Dg')
end