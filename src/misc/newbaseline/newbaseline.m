datapath = '/Users/nicthib/Documents/PhD/Data/H_92_RAW';

load(fullfile(datapath,'cm125_1_runB_H.mat'));

% Approach 1: single baseline
TC = H.red(10,:);
% add trend
TC = TC+linspace(0,400,numel(TC));
rotf = m.rotf;
blepoch = 900:1000;
t = (1:numel(rotf))/20;
close all

figure;
subplot(311)
plot(t,TC);
hold on
plot(t(blepoch),TC(blepoch))
legend('Raw Timecourse','Baseline Epoch')
title('Raw Red LED')
ylabel('counts')

subplot(312)
plot(t,rotf)
title('Movement')
subplot(313)
BL = mean(TC(blepoch));
plot([0 t(end)],[0 0],'k--'); hold on
plot(t,100*(1-TC/BL),'Color',cmap(1,:))
title('% change')
xlabel('time (sec)')
ylabel('% change')

%% Approach 2: dynamic baseline
TC = H.red(10,:);
TC = TC+linspace(0,2000,numel(TC));

blW = 100; % Baseline width
blB = 200; % Baseline buffer
blS = 1000; % Baseline skip
blT = []; blV = []; blIDX = [];
for i = blB+1:blS:numel(rotf)-blW-blB
    idx = i-blB:i+blW+blB;
    if max(rotf(idx))<1
        blT(end+1) = idx(1);
        blV(end+1) = mean(TC(idx));
        blIDX(:,end+1) = i:i+blW;
    end
end

%%
close all
figure;
subplot(311)
plot(t,TC);
hold on
plot(t(blIDX),TC(blIDX),'red')
legend('Raw Timecourse','Baseline Epochs')

subplot(312)
plot(t,rotf)
title('Movement')

subplot(313)
BL = interp1(blT,blV,1:numel(rotf),'linear','extrap');
plot([0 t(end)],[0 0],'k--'); hold on
plot(t,100*(1-TC./BL),'Color',cmap(1,:))

title('% change')
xlabel('time (sec)')
ylabel('% change')


