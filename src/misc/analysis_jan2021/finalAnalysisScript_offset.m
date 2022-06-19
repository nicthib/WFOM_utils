clear;clc;

windowlength = 60; % 3s
windowstep = 30;
numstates = 9;

mouse = 'cm128.mat';

probabilityDIR = strcat('/home/bahar/Documents/Scripts/State-to-Behavior/Bahar/2021-01-16/results/Probability/',num2str(windowlength/20),'s/',...
    num2str(numstates),'/');

load(strcat(probabilityDIR,mouse))

% maxDuration = 0;
% for i = 1:length(runningBlocks)
%     if runningBlocks(i).duration > maxDuration
%         maxDuration = runningBlocks(i).duration;
%         run = runningBlocks(i).run;
%         day = runningBlocks(i).day;
%     end
%
% end

%% sort the running blocks based on the duration

runningBlocksTable = struct2table(runningBlocks);
runningBlocksTableSorted = sortrows(runningBlocksTable, 'duration', 'descend');

% remove the runs with duration less than 5 (s)
runningBlocksTableSorted(runningBlocksTableSorted.duration < 100,:) = [];


%%
mytemp = runningBlocksTableSorted.rotf_ave(1);
maxLength = length(mytemp{1});

%% extract all rotf_ave
rotf_ave = [];
for i = 1: height(runningBlocksTableSorted)
    mytemp = runningBlocksTableSorted.rotf_ave(i);
    rotf = mytemp{1};
    if length(rotf) < maxLength
        rotf = [repelem(0,maxLength - length(rotf))';rotf];
    end
    rotf_ave = [rotf_ave rotf];
    clear mytemp rotf
end
imagesc(rotf_ave')
set(gcf, 'Position', [615   487   748   594])
colorbar
%set(gca,'XTick',(1:1:maxLength))
set(gca,'XTick',(1:2:maxLength), 'XTickLabel',(-7:2:maxLength-7),'FontSize', 9)
xlabel('Time Window','FontSize',14,'FontWeight','bold')
ylabel('Running Blocks','FontSize',14,'FontWeight','bold')
title('Running Signals Sorted based on Duration','FontSize',16)
colorbar
%caxis([0 1])
colormap jet

%% extract all states
saveDIR = '/home/bahar/Documents/Scripts/State-to-Behavior/Bahar/2021-01-16/figures/offset/cm128/';
probabilityVectorAve = zeros(maxLength, numstates);
for state = 1:numstates
    probabilityVector = [];
    for i = 1: height(runningBlocksTableSorted)
        mytemp = runningBlocksTableSorted.probabilityVector(i);
        prob = mytemp{1}(:,state);
        if length(prob) < maxLength
            prob = [repelem(0,maxLength - length(prob))';prob];
        end
        probabilityVector = [probabilityVector prob];
        clear mytemp prob
    end
    probabilityVectorAve(:,state) = mean(probabilityVector,2);
    imagesc(probabilityVector')
    set(gcf, 'Position', [615   405   837   676])
    % (maxLength * windowstep * 0.05) - 7 - 1
    %set(gca,'XTick',(1:2:maxLength), 'XTickLabel',(-7:3:68))
    set(gca,'XTick',(1:2:maxLength), 'XTickLabel',(-7:2:maxLength-7), 'FontSize', 9)
    xlabel('Time Window','FontSize',14,'FontWeight','bold')
    ylabel('Running Blocks','FontSize',14,'FontWeight','bold')
    title(strcat('Probability of State',{' '}, num2str(state)'),'FontSize',16)
    colorbar
    caxis([0 1])
    colormap jet
    pause(2)
    saveas(gcf,strcat(saveDIR,'probe_state_',num2str(state),'.tif'))
end


%% extract all pupil_ave
pupil_ave = [];
for i = 1: height(runningBlocksTableSorted)
    mytemp = runningBlocksTableSorted.pupil_ave(i);
    pupil = mytemp{1};
    if length(pupil) < maxLength
        pupil = [repelem(pupil(1),maxLength - length(pupil))';pupil];
        %pupil = [repelem(0,maxLength - length(pupil))';pupil];
    end
    pupil_ave = [pupil_ave pupil];
    clear mytemp pupil
end
pupilSignalAve = mean(pupil_ave,2);
imagesc(pupil_ave')
set(gcf, 'Position', [615   487   748   594])
colorbar
%set(gca,'XTick',(1:1:maxLength))
set(gca,'XTick',(1:2:maxLength), 'XTickLabel',(-7:2:maxLength-7),'FontSize', 9)
xlabel('Time Window','FontSize',14,'FontWeight','bold')
ylabel('Running Blocks','FontSize',14,'FontWeight','bold')
title('Pupil Signals Sorted based on Duration','FontSize',16)
colorbar
%caxis([0 1])
colormap jet

%% extract all whisk_ave
whisk_ave = [];
for i = 1: height(runningBlocksTableSorted)
    mytemp = runningBlocksTableSorted.whisk_ave(i);
    whisk = mytemp{1};
    if length(whisk) < maxLength
        whisk = [repelem(whisk(1),maxLength - length(whisk))';whisk];
        %whisk = [repelem(0,maxLength - length(whisk))';whisk];
    end
    whisk_ave = [whisk_ave whisk];
    clear mytemp whisk
end
whiskSignalAve = mean(whisk_ave,2);
imagesc(whisk_ave')
set(gcf, 'Position', [615   487   748   594])
colorbar
%set(gca,'XTick',(1:1:maxLength))
set(gca,'XTick',(1:2:maxLength), 'XTickLabel',(-7:2:maxLength-7),'FontSize', 9)
xlabel('Time Window','FontSize',14,'FontWeight','bold')
ylabel('Running Blocks','FontSize',14,'FontWeight','bold')
title('whisk Signals Sorted based on Duration','FontSize',16)
colorbar
%caxis([0 1])
colormap jet

%% plot the ave probability 
h = plot(probabilityVectorAve, 'LineWidth',1.5);
set(gcf,'Position',[615   501   728   580]);
set(h, {'color'}, {[1 0.5 0.5];[1 0 1];[0 1 1];[1 0 0];[0 1 0];[0 0 1];[0 0 0];[0.75 0.5 0.5];[0.5 0.5 1]});
Legend=cell(numstates,1);
 for iter=1:numstates
   Legend{iter}=strcat('State', num2str(iter));
 end
 legend(Legend, 'Location','best')
set(gca,'XTick',(1:2:maxLength), 'XTickLabel',(-7:2:maxLength-7),'FontSize', 9)
legend off
xline(36,'--')
xline(37,'--')
xline(38,'--')
xlim([25,42])
%ylim([0 1])
ylim([-1 1.1])
xlabel('Time Window','FontSize',14,'FontWeight','bold')
ylabel('Running Blocks','FontSize',14,'FontWeight','bold')
title('Average Probabilities Across Running Blocks','FontSize',16)
hold on
plot(pupilSignalAve/max(abs(pupilSignalAve)), 'LineWidth',1.5)
plot(whiskSignalAve/max(abs(whiskSignalAve)), 'LineWidth',1.5)





