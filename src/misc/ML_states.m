%% Initialize
clear
Hdir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; % The base folder we are pulling data from
savedir = '/local_mount/space/dingus/1/RS_analysis/analysis_march2021/d_files';
fundir = '/local_mount/space/dingus/1/RS_analysis/code';
addpath(genpath(fundir));
load('/local_mount/space/dingus/1/RS_analysis/Draft/params.mat')
% Set options
for i = 1:5
    runnames_singlemouse{i} = runnames_RS(cellfun(@(s) ~isempty(regexp(s,['cm12' mat2str(i+3)])),runnames_RS)); % Parse out runs that are 'mouse'
end
runnames_label = [runnames_singlemouse{1}(1:3:15);...
                  runnames_singlemouse{2}(1:3:15);...
                  runnames_singlemouse{3}(1:3:15);...
                  runnames_singlemouse{4}(1:3:15);
                  runnames_singlemouse{5}(1:3:15)];
              
clear results; opts = [];

%% Main params
opts.k = 7; 
opts.ww = 201; 

opts.nreps = 100; opts.sig = 0;
opts.skp = (opts.ww-1)/2;
opts.Hdir = Hdir; opts.variable = 'jrgeco';
opts.dokmeans = 1; opts.labelstates = 0;
savepath = ['/local_mount/space/dingus/1/RS_analysis/analysis_march2021/d_files/k' mat2str(opts.k) '_t' round(mat2str((opts.ww-1)/20))];
mkdir(savepath)

[~,cnt_label] = getcorrstates_v2(runnames_label,opts);   
opts.labelstates = 1;
[~,I] = sort(mean(cnt_label,2)); opts.cnt = cnt_label(I,:);

for n = 1:5
    [results{n}] = getcorrstates_v2(runnames_singlemouse{n},opts);
end

%% Combine, save
for i = 1:numel(results)
    for j = 1:numel(results{i})
        if isfield(behavior{i}{j},'whisk') && isfield(behavior{i}{j},'pupil')
            D = results{i}{j}.d;
            P = behavior{i}{j}.pupil;
            Pscore = behavior{i}{j}.p_score;
            W = behavior{i}{j}.whisk;
            Wscore = behavior{i}{j}.w_score;
            if Pscore > .95
                save(fullfile(savepath,[runnames_singlemouse{i}{j} '.mat']),'D','P','W')
            end
        end
    end
end

%save(fullfile(savepath,['cm12' mat2str(n+3) '_states.mat']),'c')



%% OLD
input_full = []; output_full = []; opts.netsize = []; opts.nepochs = 1000;
IW = [];
perf = [];
for n = 1:4
    input_full{n} = []; output_full{n} = [];
    for i = 1:20
        load(runnames_singlemouse{n}{i},'m');
        if numel(m.pupil) > 10 && numel(m.whisk) > 10
            input = d{n}{i};
            output = double([m.rotf;m.whisk;m.pupil']');
%             for j = 1:3
%                 output(:,j) = smooth(padarray(output(1:end-(opts.ww-1)/2,j),[(opts.ww-1)/2 0],'pre'),opts.ww);
%                 output(:,j) = padarray(output(1:end-(opts.ww-1)/2,j),[(opts.ww-1)/2 0],'pre');
%                 output(:,j) = smooth(output(:,j),opts.ww);
%             end
%             input_full{n} = [input_full{n}; input];
%             output_full{n} = [output_full{n}; output];
            
            for b = 1:3
                [net,~] = WFOM_MLtrain(input',output(:,b)',opts);
                IW(:,end+1,b) = net.IW{1};
                guess(b,:) = net(d{n}{i}');
            end
            output(isnan(output)) = min(output(:));
            guess(isnan(guess)) = min(guess(:));
            perf(end+1,:) = diag(corr(guess',output));
        end
    end
end

%%
net = [];
for n = 1:4
    for b = 1:3
        opts.netsize = [];
        opts.nepochs = 1000;
        [net{n}{b},~] = WFOM_MLtrain(input_full{n}',output_full{n}(:,b)',opts);
    end
end

%%
figure
for n = 1:4
    perf = [];
    for i = 1:20
        load(runnames_singlemouse{n}{i},'m','B');
        %m.whisk = 
        if numel(m.pupil) > 10 && numel(m.whisk) > 10
            output = double([m.rotf;m.whisk;m.pupil']');
            output = padarray(output(1:end-(opts.ww-1)/2,:),[(opts.ww-1)/2 0],'pre');
            for b = 1:3
                guess(b,:) = net(d{n}{i}');
            end
            output(isnan(output)) = min(output(:));
            guess(isnan(guess)) = min(guess(:));
            perf(end+1,:) = diag(corr(guess',output));
            
            close all
            titles = {'Movement','Whisking','Pupil'};
            tiledlayout(3,1);
            for b = 1:3
                ax(b) = nexttile;
                plot(output(:,b))
                hold on
                plot(guess(b,:))
                title(titles{b})
            end
            linkaxes(ax,'x')
            
        end
    end
    hold on; errorbar([1:3]+n*.05-.1,mean(perf,1),CI95(perf))
end
legend

%%
n = 1; p = 12;
load(runnames_singlemouse{n}{p},'m','B');
output = double([m.rotf;m.whisk;m.pupil']');
output = padarray(output(1:end-(opts.ww-1)/2,:),[(opts.ww-1)/2 0],'pre');
for b = 1:3
    guess(b,:) = net{n}{b}(d{n}{p}');
end
close all
titles = {'Movement','Whisking','Pupil'};
tiledlayout(3,1);
for i = 1:3
    ax(i) = nexttile;
    plot(output(:,i))
    hold on
    plot(guess(i,:))
    title(titles{i})
end
linkaxes(ax,'x')

c