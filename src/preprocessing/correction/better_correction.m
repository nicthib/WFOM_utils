%% Determine objective function(s) using gradient of converted jrgeco data
clear
addpath('/local_mount/space/dingus/1/RS_analysis/code/utils')
load('/local_mount/space/dingus/1/RS_analysis/preprocessing/params.mat', 'sessions_for_analysis')
savepath = '/local_mount/space/dingus/1/RS_analysis/preprocessing/correction';
for i = 1:numel(sessions_for_analysis)
    close all
    mouse = sessions_for_analysis{i};
    fulldir = findmousefolder(mouse,'runC','1');
    
    [m,data] = LoadData_v2(fulldir,m);
    m.aux = []; m.DAQ = []; % free up some space
    % Get some movement and some RS epoch
    m.ep1 = baselinefromrot(m.rotf,500); 
    m.ep2 = find(m.rotf > 2); m.ep2 = m.ep2(1:min([numel(m.ep2),500]));
    % combine
    fn = fieldnames(data);
    for f = 1:numel(fn)
       data.(fn{f}) = data.(fn{f})(:,:,[m.ep1 m.ep2]); 
    end
    [m.Q,m.bDr,m.bDg] = DrDg_brute_search(data,m);
    save(fullfile(savepath,[mouse '_results.mat']),'m')
    disp(['Best Dr = ', mat2str(m.bDr),' Best Dg = ', mat2str(m.bDg)])
end

%% Show results
% Gather results
clear
loadpath = '/local_mount/space/dingus/1/RS_analysis/analysis_dec2020/correction_analysis';
load('/local_mount/space/dingus/1/RS_analysis/analysis_dec2020/params.mat','sessions_for_analysis')
for i = 1:numel(sessions_for_analysis)
    load(fullfile(loadpath,[sessions_for_analysis{i} '_results.mat']))
    Dr_all(i) = m.bDr;
    Dg_all(i) = m.bDg;
    Q(:,:,i) = m.Q;
end
%Q(:,:,[5 6 7 11 12 21]) = [];
Qm = mean(Q,3);
Dri = m.Drt; Dgi = m.Dgt;
Qmin = min(Qm(:));
[R,C] = find(Qm==Qmin);
bDr = Dri(R); bDg = Dgi(C);

% Show mean objective function
close all
figure(1)
[X,Y] = meshgrid(Dri,Dgi);
scatter3(bDr,bDg,min(Qm(:))+.1,25,'filled','MarkerFaceColor','r')
hold on
surf(X,Y,Qm')
shading interp
view(2)
xlabel('X_r')
ylabel('X_g')
zlabel('Objective function')
colormap jet
title('max(G(F(t,X_r,X_g)))')
axis image; colorbar
text(bDr,bDg,.2,['  ' mat2str(bDr) ', ' mat2str(bDg)],'Color','w')

% Show session objective functions
figure(2)
for i = 1:size(Q,3)
    subplot(6,4,i)
    Qtmp = Q(:,:,i);
    [~,idx] = min(Qtmp(:));
    [R,C] = ind2sub(size(Qtmp),idx);
    bDrtmp = Dri(R);
    bDgtmp = Dgi(C);
    scatter3(bDrtmp,bDgtmp,1,10,'filled','MarkerFaceColor','r')
    hold on
    surf(X,Y,Q(:,:,i)')
    shading interp
    colormap jet
    axis image
    view(2); 
    if i == 1
        xlabel('Dr'); ylabel('Dg')
    else
        axis off
    end
    title(strrep(sessions_for_analysis{i},'_',' '))
end
clear bDrtmp bDgtmp Qtmp R C
saveas(1,'Mean_obj_function.png')
saveas(2,'Indiv_obj_function.png')

%% Older method
%m.stim_data = meanstim(m,data);
%stim_mean = mean(m.stim_data.lime(:,:,100:130),3)./mean(m.stim_data.lime(:,:,1:70),3);
%imagesc(stim_mean); axis image; axis off; colormap gray; title('Pick response region')
%m.whisk_roi = double(roipoly); m.whisk_roi(m.whisk_roi==0) = NaN; close all

%fun = @(D) test_convert(data,D,m,'RS');
%options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',50);
%[D(i,:),Q(i),~,~] = fminsearch(fun,[.5,.5],options);
%figure(1); sgtitle(strrep(mouse,'_',' '));
%saveas(1,[mouse '.fig'])
%saveas(2,[mouse '_loss.fig'])