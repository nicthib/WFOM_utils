%% Initialize
clearvars -EXCEPT pths
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))
savedir = '/local_mount/space/chronos/2/RS_SVD_2020';
load('/local_mount/space/dingus/1/RS_analysis/cm_analysis/bkg.mat');

%% SVD 
runs = 'BCDEFGHIJ';
days = '12345678';
mice = {'125','126','127','128'};
for a = 1:numel(mice)
    for b = 1:numel(days)
        for c = 1:numel(runs)
            tic
            m = makem;
            m.mouse = ['cm' mice{a} '_' days(b)]; m.dsf = 2; m.nrot = 2; m.PCAcomps = 500;
            m.run = ['run' runs(c)]; m.stim = '1'; m.ccddir = findmousefolder(m.mouse);
            m.outputs = 'P'; m.BW = ones(512/m.dsf,512/m.dsf);
            m.bkgsub = bkg;
            fulldir = findmousefolder(m.mouse,m.run,m.stim);
            if ~isempty(fulldir)
                try
                    [m,data] = LoadData_v21(fulldir,m);
                    C = data.C; S = data.S; expl = data.expl;
                    save(fullfile(savedir,[m.mouse '_' m.run]),'-v7.3','C','S','expl','m')
                    disp(['Done, took ' mat2str(round(toc)) ' sec'])
                catch
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% %% Reconstruction and creation of H
% clear
% load('/local_mount/space/dingus/1/RS_analysis/H_new/IDX.mat')
% IDXs = IDX; clear IDX
% mice = {'cm125','cm126','cm127','cm128'};
% runs = 'BCDEFGHIJ';
% days = '12345678';
% savedir = '/local_mount/space/dingus/1/RS_analysis/H_18';
% cd(savedir)
% for a = 1:numel(mice)
%     for b = 1:numel(days)
%         for c = 1:numel(runs)
%             try
%                 H = []; Hstd = [];
%                 filename = [mice{a} '_' days(b) '_run' runs(c) '.mat'];
%                 savefilename = [mice{a} '_' days(b) '_run' runs(c) '_H.mat'];
%                 load(filename);
%                 rotf = getrotf(m);
%                 m.baseline = baselinefromrot(rotf,100,1);
%                 disp(['Converting ' filename])
%                 for d = 1:m.nLEDs
%                     data.(m.LEDs{d}) = convertSVD(C.(m.LEDs{d}),S.(m.LEDs{d}),1:500,m.nanidx,256);
%                 end
%                 data.jrgeco = data.lime./((data.red.^m.Dr).*(data.green.^m.Dg));
%                 [data.chbo,data.chbr,~] = convert_mariel_MIAO(data.green,data.red,'g','r',m.baseline,m.greenfilter);
%                 bgGG = mean(data.jrgeco(:,:,m.baseline),3);
%                 data.jrgeco = data.jrgeco./repmat(bgGG,[1 1 size(data.jrgeco,3)])-1;
%                 % reshape to speed up H
%                 data.jrgeco = reshape(data.jrgeco,[256^2,size(data.jrgeco,3)]);
%                 data.chbo = reshape(data.chbo,[256^2,size(data.chbo,3)]);
%                 data.chbr = reshape(data.chbr,[256^2,size(data.chbr,3)]);
%                 for d = 1:numel(IDXs)
%                     IDX = IDXs{d};
%                     if mice{a} ~= 'cm125'
%                         IDX = imwarp_fr(IDX,TF_a.(['cm125to' mice{a}]));
%                     end
%                     if days(b) ~= '1'
%                        IDX = imwarp_fr(IDX,TF_w.([mice{a} '_1_' days(b)]));
%                     end
%                     % GetHfromKmeans
%                     for e = 1:max(IDX(:))
%                         tmp = double(IDX == e);
%                         tmp(tmp == 0) = NaN; idx = ~isnan(reshape(tmp,[1 256^2]));
%                         H.jrgeco{d}(e,:) = squeeze(nanmean(data.jrgeco(idx,:),1));
%                         H.chbo{d}(e,:) = squeeze(nanmean(data.chbo(idx,:),1));
%                         H.chbr{d}(e,:) = squeeze(nanmean(data.chbr(idx,:),1));
%                         Hstd.jrgeco{d}(e,:) = squeeze(nanstd(data.jrgeco(idx,:),1));
%                         Hstd.chbo{d}(e,:) = squeeze(nanstd(data.chbo(idx,:),[],1));
%                         Hstd.chbr{d}(e,:) = squeeze(nanstd(data.chbr(idx,:),[],1));
%                     end
%                 end
%                 [st,wz,pz,rz,whisk,pupil,rotf] = getbehaviorz(savefilename);
%                 m.DAQ = []; m.aux = [];
%                 save(fullfile(savedir,savefilename),'rotf','pupil','whisk','rz','wz','pz','H','Hstd','m','-v7.3')
%             catch me
%                 disp(['failed for ' mice{a} '_' days(b) '_run' runs(c) ', error is ' me.message])
%             end
%         end
%     end
% end
% 
% %% Optimal state analysis
% %% Initialize
% clear
% root = '/local_mount/space/dingus/1/RS_analysis/H_new';
% vars = {'H'};
% allrunnames = getallexps(root,vars)'; % Gets all runs in the H_nf folder
% 
% %% Set options
% clearvars -EXCEPT allrunnames;
% d_mean = [];
% load('TF.mat','cm125')
% mice = {'cm125','cm126','cm128'};
% for m = 1:numel(mice)
%     mouse = mice{m}; runs = 'BCD';
%     opts.mapidx = 1;  opts.nreps = 500; opts.sig = 6;
%     opts.ww = 41; opts.skipfactor = (opts.ww-1)/2;
%     o = cm125.IDXparts{opts.mapidx}; % For visual division of corr maps
%     divs = [numel(o.L), numel([o.L o.C])]; % For visual division of corr maps
%     runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames))'; % Parse out runs that end with B,C,D
%     runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames))'; % Parse out runs that are 'mouse'
%     runnames_test = runnames(1:3:end); % Get subset
%     
%     d_mean.(mouse)=[];
%     d_std.(mouse)=[];
%     for k = 1:20
%         opts.k = k;
%         opts.dokmeans = 1; opts.labelstates = 1;
%         [~,~,~,d] = getcorrstates(runnames_test,opts);
%         dfull = cell2mat(d');
%         d_mean.(mouse)(k) = mean(dfull(:));
%         d_std.(mouse)(k) = std(dfull(:));
%     end
% end
% 
% %% 
% 
% for i=2:18
%     LM1 = fitlm(2:i,d_mean.cm125(2:i));
%     LM2 = fitlm(i+1:20,d_mean.cm125(i+1:20));
%     R(i) = (LM1.Rsquared.Ordinary+LM2.Rsquared.Ordinary)/2;
%     
% end

