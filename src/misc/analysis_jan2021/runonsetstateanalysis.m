%% Initialize
clear
root = '/local_mount/space/dingus/1/RS_analysis/H_new';
vars = {'H','rotf','whisk','pupil'};
allrunnames = getallexps(root,vars)'; % Gets all runs in the H_nf folder

%%
clearvars -EXCEPT allrunnames;
load('IDX.mat')

%% Set options and run names
opts.mapidx = 2; opts.k = 12; opts.nreps = 100; opts.sig = 1;
opts.ww = 33; opts.skipfactor = 20;
mice = {'cm125','cm126','cm128'};
o = IDXparts{opts.mapidx}; % For visual division of corr maps
divs = [numel(o.L), numel([o.L o.C])]; % For visual division of corr maps
for m = 1:numel(mice)
    mouse = mice{m};
    runnames.label = getrunnames(mouse,'run[BCD]',allrunnames);
    runnames.RS.(mouse) = getrunnames(mouse,'run[BCD]',allrunnames);
    tmp = getrunnames(mouse,'run[EFG]',allrunnames)';
    runnames.ST.(mouse) = [];
    for j = 1:numel(tmp)
        load(strrep(tmp{j},'_H',''),'m')
        if isfield(m,'stimframes')
            runnames.ST.(mouse){end+1,1} = tmp{j};
        end
    end
    win = getwin(opts.ww,opts.sig);
    twr = (numel(win)-1)/2;
    
    % Get states
    opts.dokmeans = 1; opts.labelstates = 0;
    [~,c] = getcorrstates(runnames.label(1:2:end),opts);
    [~,I] = sort(sum(c,2),1,'ascend');
    c = c(I,:);
    
    % Label states
    opts.dokmeans = 0; opts.labelstates = 1; opts.cnt = c;
    [~,~,st.RS.(mouse)] = getcorrstates(runnames.RS.(mouse),opts);
    [~,~,st.ST.(mouse)] = getcorrstates(runnames.ST.(mouse),opts);
    disp(['DONE ' mouse])
end


%%
st_on = []; rotf_on = []; pupil_on = []; whisk_on = []; 
stimruns = 0; mouse = 'cm125';
if stimruns == 0
    runnamestotest = runnames.RS.(mouse);
    sttotest = st.RS.(mouse);
else
    runnamestotest = runnames.ST.(mouse);
    sttotest = st.ST.(mouse);
end
rslength = 200;
for i = 1:numel(runnamestotest)
    load(runnamestotest{i},'pupil','rotf','whisk')
    if stimruns
        if isfield(m,'stimframes')
            onframes = m.stimframes;
        else
            onframes = [];
        end
        for j = 2:numel(onframes)
            idx=onframes(j)-rslength:onframes(j)+rslength;
            st_on(end+1,:,:) = sttotest{i}(idx,:);
            rotf_on(end+1,:) = rotf(idx);
            pupil_on(end+1,:) = zscore(pupil(idx));
            whisk_on(end+1,:) = whisk(idx);
        end
    else
        [onframes,offframes] = getrunonsetframes(runnamestotest{i},0);
        for j = 2:numel(onframes)
            if onframes(j)-offframes(j-1) >= rslength
                idx=onframes(j)-rslength:onframes(j)+rslength;
                if min(idx) > 0 && max(idx) < 11980
                    st_on(end+1,:,:) = sttotest{i}(idx,:);
                    rotf_on(end+1,:) = rotf(idx);
                    pupil_on(end+1,:) = zscore(pupil(idx));
                    whisk_on(end+1,:) = whisk(idx);
                end
            end
        end
    end
end

%
st_dist = []; st_tmp = [];
for i = 1:opts.k
    st_tmp = squeeze(st_on(:,:,i));
    st_dist(i,:) = mean(st_tmp);
    st_std(i,:) = std(st_tmp);
end

a=figure;
leg = [];
for s = 1:2
    spstates = [1:size(st_dist,1)/2]+6*(s-1);
    for l = 1:opts.k/2
       leg{l} = ['State ' mat2str(spstates(l))]; 
    end
    subplot(5,1,[1:2]+(s*2-2))
    hold on
    if stimruns
        txt = 'Stimulus';
    else
        txt = 'Spontaneous';
    end
    if s == 1
        title([mouse ' ' txt])
    end
    cmap = jet(size(st_dist,1))*.8;
    t=linspace(-rslength/20,rslength/20,size(st_dist,2));
    t=t+twr/20;
    if stimruns
        %stimframes are off by 5 frames
        t = t-5/20;
    end
    for i = spstates
        tmp = st_dist(i,:);
        loc = find(tmp==max(tmp));
        plot(t,tmp,'Color',cmap(i,:),'Linewidth',2);
        
        text(t(loc(1)),max(tmp)+1,num2str(i),'Color',cmap(i,:))
    end
    stdshade(st_std,.2,cmap(i,:),t)
    legend(leg,'Location','EastOutside')
    plot([0 0],[0 1000],'Color','k')
    ylabel('Avg. Distance from state centroid')
    xlim([-10 10])
    grid on
end
p_plot = mean(pupil_on);
r_plot = mean(rotf_on);
subplot(5,1,5)
plot(t-twr/20,p_plot,'Linewidth',2)
hold on
plot(t-twr/20,r_plot,'Linewidth',2)
xlim([-10 10])
ylim([-2 2])
legend('pupil','rotary','Location','EastOutside')
xlabel(['Time since ' txt ' (sec)'])

%% Other visualize
b=figure;
hold on
ctmp = c;
r = sqrt(size(ctmp,2));
k = size(ctmp,1);
for j = 1:k
    if j <= 6
        spn = j;
    else
        spn = j+6;
    end
    subplot(4,6,spn)
    tmp = divcorrmap(reshape(ctmp(j,:),[r r]),divs);
    im = imagesc(tmp);
    im.AlphaData = ~isnan(tmp);
    caxis([0 1]);
    axis image; axis off
    title(['State #' mat2str(j)])
    colormap(jet)
    subplot(4,6,spn+6)
    tmp = corrIDXmap(IDX{opts.mapidx},reshape(ctmp(j,:),[r r]));
    tmp(tmp == 0) = NaN;
    im = imagesc(tmp);
    im.AlphaData = ~isnan(tmp);
    axis image; caxis([0 1]); axis off
end

%%
figure
subplot(211); hold on
title('Pupil')
plot(t-(twr)/20,p_plot1,'Color',[0 0 .75])
plot(t-(twr)/20,p_plot2,'Color',[.5 .5 1])
legend('Sleepy stim','Waking stim')
xlim([-10 10])
ylim([-2 1])
subplot(212); hold on
title('Movement')
plot(t-(twr)/20,r_plot1,'Color',[.75 0 0])
hold on
plot(t-(twr)/20,r_plot2,'Color',[1 .5 .5])
legend('Sleepy stim','Waking stim')
xlim([-10 10])
ylim([-2 1])


%% Run offset
st_on = cell(opts.k,1); rotf_on = cell(opts.k,1); pupil_on = cell(opts.k,1); 
mouse = 'cm128';
runnamestotest = runnames.RS.(mouse);
sttotest = st.RS.(mouse);
st_off = zeros(1,opts.k,1); rotf_off = []; pupil_off = [];
rslength = 200;
maxoff = 0;
for i = 1:numel(runnamestotest)
    [on,off] = getrunonsetframes(runnamestotest{i},rslength,0);
    if max(off(2:end)-on(1:end-1)) > size(st_off,3)
        maxoff = max(off(2:end)-on(1:end-1));
        st_off = padarray(st_off,[0 0 maxoff-size(st_off,3)],'post');
        rotf_off = padarray(rotf_off,[0 maxoff-size(rotf_off,2)],'post');
        pupil_off = padarray(pupil_off,[0 maxoff-size(pupil_off,2)],'post');
    end
    load(runnamestotest{i},'rotf','pupil');
    for j = 1:numel(off)-1
        if on(j+1)-off(j) > rslength
            idx=off(j):on(j+1);
            st_off(end+1,:,:) = cat(2,zeros(opts.k,maxoff-numel(idx)),sttotest{i}(idx,:)');
            rotf_off(end+1,:) = [zscore(rotf(idx))' zeros(1,maxoff-numel(idx))];
            pupil_off(end+1,:) = [zscore(pupil(idx))' zeros(1,maxoff-numel(idx))];
        end
    end
end

%%
st_off(st_off==0) = NaN;
st_dist = [];
for i = 1:opts.k
   st_tmp = squeeze(st_off(:,i,:));
   st_dist(i,:) = nanmean(st_tmp,1);
end
