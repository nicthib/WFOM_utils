%% Setup
clear;
opts.ww = 45; opts.mapidx = 1; opts.k = 10; opts.nreps = 500;
opts.skipfactor = floor(opts.ww/2);
opts.sig = 6; opts.dokmeans = 1;
runnames = getallexps('/local_mount/space/dingus/1/RS_analysis/analysis_dec2020/H',{'H'}); % Gets all runs in the H_nf folder
runnames = runnames(cellfun(@(s)~isempty(regexp(s,'run[BCD]')),runnames))'; % Parse out runs that end with B,C,D
%runnames = runnames(cellfun(@(s)~isempty(regexp(s,'125_6')),runnames))'; % Parse out runs that are named blah
bootstrap = 5; % # of runs that will be randomly selected per bootstrap

%% Get states
% for i = 1:3 % BOOTSTRAPS
%     runidx = randperm(numel(runnames),bootstrap);
%     [st{i},c{i},d{i}] = getcorrstates({runnames{runidx}},opts);
% end

% NO BOOTSTRAP
i=1;
[st{i},c{i},d{i}] = getcorrstates(runnames,opts);

r = sqrt(size(c,2)); k = opts.k;
disp('DONE!')

%% Show corr maps
if iscell(c) && numel(c) > 1
    [ord, distmin] = orderstates(c,1:k);
else
    ctmp = c; c = [];
    c{1} = ctmp;
    ord = 1:k;
end
f1 = figure;
for i = 1:size(c,2)
    for j = 1:k
        figure(f1)
        subplot(size(c,2),k,(i-1)*k+ord(i,j))
        imagesc(reshape(c{i}(j,:),[r r]))
        caxis([0 1]); 
        axis image; axis off
        %title(['S' num2str(j), ', dist = ' mat2str(round(distmin(i,j)/r,2))])
    end
end
colormap(jet)

%% Plotting windows
win1 = getwin(opts.sp,opts.sig);
win2 = getwin(opts.sp,0);
pad = (numel(win1)-numel(win2))/2;
close all
tmp = (numel(win1)-1)/2;
t = -tmp:tmp;
plot(t,[zeros(1,pad) win2 zeros(1,pad)])
hold on
plot(t,win1)
ylim([-1 2])
legend('Square window','Tapered window')
xlabel('Frame #')

%% Link orders (for long vs short)
[ordc, ~] = orderstates(c,1:k);
for i = 1:numel(c2)
    c2{i} = c2{i}(ordc(2,:),:);
end
% combine
c = {c1{:},c2{:}};

