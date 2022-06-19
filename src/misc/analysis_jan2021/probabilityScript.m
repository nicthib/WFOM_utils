clear;clc;

windowlength = 60; % 3 seconds
windowstep = 30;
numstates = 9;

mouse = 'cm128.mat';
runningBlocksDIR = '/home/bahar/Documents/Scripts/State-to-Behavior/Bahar/2021-01-16/results/runningBlocks/';
stateDIR = strcat('/home/bahar/Documents/Scripts/State-to-Behavior/Bahar/2021-01-09/results/States/',num2str(windowlength/20),'s/',...
    num2str(numstates),'/');
saveDIR = strcat('/home/bahar/Documents/Scripts/State-to-Behavior/Bahar/2021-01-16/results/Probability/',num2str(windowlength/20),'s/',...
    num2str(numstates),'/');

load(strcat(runningBlocksDIR,mouse))

% calculate correlation windows for each run
for run = 1:length(runningBlocks)
    % the name of neural states associated with the run
    statename = strcat(runningBlocks(run).mouse_name,'_',runningBlocks(run).day,'_',runningBlocks(run).run,'.mat');
    
    % load the states
    load(strcat(stateDIR,statename))
    
    % calculate correlation windows for each run
    neuralsig = runningBlocks(run).neural;
    neuralsig = neuralsig';
    rotf = runningBlocks(run).rotf;
    pupil = runningBlocks(run).pupil;
    whisk = runningBlocks(run).whisk;
    
    
    ss = size(neuralsig);
    
    lim = ss(2) - windowlength;
    endpoint = floor(ss(2)/windowstep) * windowstep;
    
    while endpoint > lim
        endpoint = endpoint - windowstep;
    end
    
    m=0;
    %mappy = zeros(ceil(ss(2)/10),ss(1),ss(1));
    
    for i = 1:windowstep:endpoint
        m=m+1; % number of windows
        
        % save average rotf, pupil, and whisk within each window
        rotf_ave(m,1) = mean(rotf(i + [0:windowlength],1));
        pupil_ave(m,1) = mean(pupil(i + [0:windowlength],1));
        whisk_ave(m,1) = mean(whisk(i + [0:windowlength],1));
        
        for j = 1:ss(1)
            for k = 1:ss(1)
                c = corrcoef(neuralsig(j,i+[0:windowlength]),neuralsig(k,i+[0:windowlength]));
                mappy(m,j,k)=c(1,2);
            end
        end
        disp(i);
    end
    
    runningBlocks(run).correlations = mappy;
    
    % Convert correlation matrices to correlation vectors
    % use the whole matrix
    mappyvec = zeros(m,ss(1)*ss(1));
    for window = 1:m
        mytemp1 = squeeze(mappy(window,:,:));
        mytemp2 = mytemp1 -diag(diag(mytemp1));
        mappyvec(window,:) = reshape(mytemp2,[1,ss(1)*ss(1)]);
    end
    
    % calculate the probability of belonging each window to the states
    probability = [];
    for window = 1:m
        V1 = mappyvec(window,:);
        V1 = V1/norm(V1);
        for k = 1:numstates
            V2 = C_sorted(k,:);
            V2 = V2/norm(V2);
            probability_vector(1,k) = (V1 * V2')^2;
            clear V2
        end
        clear V1
        probability_matrix(window,:,:) = reshape(probability_vector,sqrt(numstates),sqrt(numstates))';
        probability = [probability;probability_vector];
    end
    runningBlocks(run).probabilityMatrix = probability_matrix;
    runningBlocks(run).probabilityVector = probability;
    
    % to save the average of running, whisking, and pupil withing each window
    runningBlocks(run).rotf_ave = rotf_ave;
    runningBlocks(run).pupil_ave = pupil_ave;
    runningBlocks(run).whisk_ave = whisk_ave;
    
    clearvars -except mouse numstates runningBlocks runningBlocksDIR saveDIR stateDIR windowlength windowstep
end
save(strcat(saveDIR,mouse),'runningBlocks')