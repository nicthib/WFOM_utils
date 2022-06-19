%% v1 
% (1) motion energy PC in differnt parts of body
% (2) camera data (in different file)

%% Add these path
clear all; close all; clc;
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
addpath(genpath('/local_mount/space/dingus/1/Weihao/Resting_State/Utils'));
addpath(genpath('/local_mount/space/dingus/1/RS_analysis'))
save_root = '/local_mount/space/dingus/1/RS_analysis/webcams/nose';

%% Pool one mouse at a time
mouse = 'cm125';
days = ['1','2','3','4','5','6'];

%% 1- which camera is the face camera?? (do it automatically)
cam_id = 0;
clear config_file img
for d = 1:numel(days)
    mouseday = [mouse '_' days(d)];
    mdir = findmousefolder(mouseday);
    webdir = dir(fullfile(mdir,'webcam','run*'));
    for i_run = 1:numel(webdir)
        img(d).(webdir(i_run).name) = LoadFLIR(fullfile(webdir(i_run).folder,webdir(i_run).name),1,cam_id);
        if ~isempty(img(d).(webdir(i_run).name))
            if size((img(d).(webdir(i_run).name)),1) == 540
                config_file(d).(webdir(i_run).name) = 1;
            else
                config_file(d).(webdir(i_run).name) = 0;
            end
        end
    end 
end

%% now we have correct cam_id
clear img
for d = 1:numel(days)
    mouseday = [mouse '_' days(d)];
    mdir = findmousefolder(mouseday);
    webdir = dir(fullfile(mdir,'webcam','run*'));
    for i_run = 1:numel(webdir)
        try
        if ~isempty(config_file(d).(webdir(i_run).name))
            cam_id = config_file(d).(webdir(i_run).name);
            img(d).(webdir(i_run).name) = LoadFLIR(fullfile(webdir(i_run).folder,webdir(i_run).name),1,cam_id);
        end
        catch
            continue
        end
    end 
end

%% crop (do it manually)
areas = {'nose'};

clear cam x y xx yy;
for d = 1:numel(days)
    
    close all; 
    figure('Position',[0 0 2000 900]);
    for i = 1:numel(areas)
        % crop on runB image (same crop for the same day)
        subplot(1,numel(areas),i)
        imagesc(img(d).(webdir(2).name));colormap gray;axis image;  title(['crop  ' areas{i}])
        [y, x] = ginput(2);
        x = round(x);
        y = round(y);
        % save
        xx(d).(areas{i})=x;
        yy(d).(areas{i}) = y;
        hold on
        plot([min(y) max(y) max(y) min(y) min(y)],[min(x) min(x) max(x) max(x) min(x)],'--y','LineWidth',2);
        
    end
    if ~exist(fullfile(save_root,mouse,days(d)))
        mkdir(fullfile(save_root,mouse,days(d)))
    end
    saveas(gcf,fullfile(save_root,mouse,days(d),'webcam_crop.jpg'));
end
%% save config file

save(fullfile(save_root,mouse,'config'),'xx','yy','img','areas','config_file','days','mouse');
clear data chbt_reduced jrgeco_reduced
%% Now let's corp

runs = fieldnames(config_file);
bug_run = [];
for d = 1:numel(days)
    for i_run = 1:10 % runA to runJ
        % load webcam in
        dsf = 1; % currently only support no downsample version
        if ~isempty(config_file(d).(webdir(i_run).name))
            mouseday = [mouse '_' days(d)];
            cam_id = config_file(d).(webdir(i_run).name);
            run = runs{i_run};
            try

            %% crop now
            disp(['done loading ' mouseday,run])
            
            for i = 1:numel(areas)
                x = xx(d).(areas{i});
                y = yy(d).(areas{i});
                cam.(areas{i}) = LoadcropFLIRmovie_wx(mouseday,run,cam_id,dsf,x,y);
            end
            %% creative part
            % originally this part is for motion energy pca after cropping
            % the video. now adapt to just saving cropped video
            
            for i = 1:numel(areas)
                % convert gray image to rgb ( basically replicate 3 times
                % of gray)
                
                rgb_movie = repmat(cam.(areas{i}),[1 1 1 3]);
                
                % Write sample frames to disk
                save_name = fullfile(save_root,mouse,days(d),[run '_nose.avi']);
                v = VideoWriter(save_name);
                v.open;
                for t = 1:size(cam.(areas{i}),3) % save frame by frame
                    v.writeVideo(uint8(squeeze(rgb_movie(:,:,t,:))));
                end
                v.close;
                disp(['done saving' mouseday,run]);
            end
            
       
% % % % % %             %% Motion energy PCA
% % % % % %             clear movie
% % % % % %             npc = 100; % keep 100 comps
% % % % % %             for i = 1:numel(areas)
% % % % % %                 motion.(areas{i}) = cam.(areas{i})(:,:,2:end) - cam.(areas{i})(:,:,1:end-1);
% % % % % %                 ss_motion = size(motion.(areas{i}));
% % % % % %                 motion.(areas{i}) = reshape(motion.(areas{i}),[ss_motion(1)*ss_motion(2) ss_motion(3)]);
% % % % % %                 [motionC.(areas{i}),motionS.(areas{i}),motionExplained.(areas{i})] = fsvd(motion.(areas{i}),npc,2,true);
% % % % % %             end
% % % % % %             disp(['done me pca ' mouseday,run])
% % % % % %             %% save all
% % % % % %             %save(fullfile(save_root,mouse,days(d),[run '_cam_pc']),'motionC','motionS','motionExplained','areas','cam')
% % % % % %             save(fullfile(save_root,mouse,days(d),[run '_cam_pc']),'motionC','motionS','motionExplained','areas')
            catch ME
                bug_run(d).(webdir(i_run).name)=1;
                disp('bug')
            end
        end
    end
end
clear movie
%%
% % % % tmp = dir(fullfile(webdir,['*_cam0.jpg']));
% % % % total_frames = numel(tmp);
% % % % time_stamp = 1:3:total_frames;
% % % clear tmp
% % % cam_id = 0;
% % % % pre-allocate
% % % tmp = LoadFLIR(webdir,1,cam_id);
% % % movie = zeros(size(tmp,1)/dsf,size(tmp,2)/dsf,numel(time_stamp));
% % % 
% % % 
% % % movie(:,:,t) = imresize(LoadFLIR(webdir,real_ind,cam_id),[size(tmp,1)/dsf,size(tmp,2)/dsf]);
% % % 
% % % 
% % % %% Load webcam
% % % tic
% % % cam_id = 0; dsf = 5;
% % % mouseday = 'cm128_3'; run = 'runC';
% % % movie1 = LoadFLIRmovie_wx(mouseday,run,cam_id,dsf);
% % % toc

