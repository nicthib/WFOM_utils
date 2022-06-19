function [auROC] = auroc(varargin)
%
%   Inputs:     If no input is given, a tutorial mode is launched
%               In tutorial mode, three graphs are generated
%               showing the extreme values of the auROC that are possible
%               for a fictitous dataset
%   
%               If any inputs are given, then three must be given and they should be
%  
%               (1) the raw data arranged as a matrix where each row corresponds to one set of observations (a time series, for example),
%               (2) an array containing the indices for the first group of data points (can be discontiguous),
%               (3) an array containing the indices for the first group of data points (can also be discontiguous, different length),
%               (4) Optional: The max for the classifier (if you believe
%               there is a maximum value for your data, this will produce
%               more consistency - where larger deviations from baseline
%               periods will produce larger auROCs)

%               in that order.
%   
%   Output:     An array of auROC values, indexed by the observations in the raw data

if nargin>0

    raw = varargin{1};
    zone1 = varargin{2};
    zone2 = varargin{3}; 

    if nargin>3
        highestpoint = varargin{4}; 
    end
    
    c_before = raw(:,zone1)';
    c_after = raw(:,zone2)';
    
    if ~exist('highestpoint'); highestpoint = max([c_before(:);c_after(:)]); end
    lowestpoint = min([c_before(:);c_after(:)]);
    increment = (highestpoint - lowestpoint)/100; 
    auROC = [];

    for j = 1:size(raw,1)

        combinedata = [c_before(:,j);c_after(:,j)];
        codes = [zeros(numel(c_before(:,j)),1);ones(numel(c_after(:,j)),1)]; 

        %loop through thresholds
        for i = 1:100
            threshold = i*increment + lowestpoint;
            result = gt(combinedata,threshold);
            %calculate true positives and false positives
            tp(i)=sum((result==1)&(codes==1))/numel(c_after(:,j));
            fp(i)=sum((result==1)&(codes==0))/numel(c_before(:,j));
        end

        % Calculate the Riemann sum underneath the curve (AUC)
        [a,b] = sort(fp); 
        auROC(j) = trapz(sort(fp),tp(b));

    end

else
    
    % HIGH AUC Figure
    
    figure('color','w'); 
    subplot(3,6,[1:3]); 
    plot([linspace(1,200,200)],[rand(1,200)],'r'); hold on;
    plot(300+linspace(1,100,100),[rand(1,100)],'k');
    plot([200+linspace(1,100,100)],[1+rand(1,100)],'b');
    title(sprintf('User can optionally specify the\nperiods of interest (red and blue)'))
    set(gca,'TickDir','out'); box off;
    subplot(3,6,[4:6]); plot(.1*(0.5-rand(1,200)+zeros(1,200)),[rand(1,100),rand(1,100)],'r.'); 
    hold on;  plot(.1*(0.5-rand(1,100))+.2*ones(1,100),[1+rand(1,100)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;
    title(sprintf('Best case scenario is perfect\n separation between red and blue'))
    widebox = [0,0,1.2,.65]*1000;
    set(gcf,'position',widebox)

    subplot(3,6,7); 
    threshold = 0.2;
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))+zeros(1,numel(data_ss))),data_ss,'r.'); 
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))+zeros(1,numel(data_ss2))),data_ss2,'b.'); 
    hold on;  plot(.1*(0.5-rand(1,100))+.2*ones(1,100),[1+rand(1,100)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,8);
    threshold = 0.6; 
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))+zeros(1,numel(data_ss))),data_ss,'r.');
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))+zeros(1,numel(data_ss2))),data_ss2,'b.'); 
    hold on;  plot(.1*(0.5-rand(1,100))+.2*ones(1,100),[1+rand(1,100)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,9);
    threshold = 1; 
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))+zeros(1,numel(data_ss))),data_ss,'r.'); 
    title('The perfect threshold');
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))+zeros(1,numel(data_ss2))),data_ss2,'b.'); 
    hold on;  plot(.1*(0.5-rand(1,100))+.2*ones(1,100),[1+rand(1,100)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,10);
    threshold = 1.2; 
    data1 = 1+[rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,200)+zeros(1,200)),[rand(1,100),rand(1,100)],'r.'); 
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss)))+.2*ones(1,numel(data_ss)),data_ss,'r.'); 
    hold on;  plot(.1*(0.5-rand(1,numel(data_ss2)))+.2*ones(1,numel(data_ss2)),data_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,11);
    threshold = 1.6; 
    data1 = 1+[rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,200)+zeros(1,200)),[rand(1,100),rand(1,100)],'r.'); 
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss)))+.2*ones(1,numel(data_ss)),data_ss,'r.'); 
    hold on;  plot(.1*(0.5-rand(1,numel(data_ss2)))+.2*ones(1,numel(data_ss2)),data_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,12);
    threshold = 1.9; 
    data1 = 1+[rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,200)+zeros(1,200)),[rand(1,100),rand(1,100)],'r.'); 
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss)))+.2*ones(1,numel(data_ss)),data_ss,'r.'); 
    hold on;  plot(.1*(0.5-rand(1,numel(data_ss2)))+.2*ones(1,numel(data_ss2)),data_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    tbrate = [1,1,1,0.9,0.4,0.1];
    fbrate = [0.9,0.4,0,0,0,0];

    colors = hsv(6);
    
    subplot(3,6,13); 
    patch([0,1,1,0],[0,0,1,1],'r','edgecolor','none','facealpha',0.2); hold on; 
    for i = 1:6
        plot(fbrate(i),tbrate(i),'.','color',colors(i,:),'markersize',18);
    end
    text(0.3,0.3,sprintf('auROC=1'))
    xlabel('False blues'); ylabel('True blues'); box off;
    set(gca,'YTick',[0:.25:1]);
    ylim([0,1.5]);

    for i = 13:18
        subplot(3,6,i);
        ylim([0,1.5]);
        text(0.1,1.45,sprintf('True blues: %1.2f',tbrate(i-12)),'color',colors(i-12,:)); 
        text(0.1,1.3,sprintf('False blues: %1.2f',fbrate(i-12)),'color',colors(i-12,:)); 
        if i~=13; set(gca,'YTick',[],'XTick',[],'YColor','w','Xcolor','w'); end
    end

    
    % 
    % LOW AUC Figure
    figure('color','w'); 
    subplot(3,6,[1:3]); 
    plot([linspace(1,200,200)],[1+rand(1,200)],'r'); hold on;
    plot(300+linspace(1,100,100),[1+rand(1,100)],'k');
    plot([200+linspace(1,100,100)],[rand(1,100)],'b');
    title(sprintf('User can optionally specify the\nperiods of interest (red and blue)'))
    set(gca,'TickDir','out'); box off;
    subplot(3,6,[4:6]); plot(.1*(0.5-rand(1,200)+zeros(1,200)),[1+rand(1,200)],'r.'); 
    hold on;  plot(.1*(0.5-rand(1,100))+.2*ones(1,100),[rand(1,100)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;
    title(sprintf('Best case scenario is perfect\n separation between red and blue'))
    widebox = [0,0,1.2,.65]*1000;
    set(gcf,'position',widebox)
    
    subplot(3,6,7); 
    threshold = 0.2;
    data1 = [rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))+2*ones(1,numel(data_ss))),data_ss,'r.'); 
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))+2*ones(1,numel(data_ss2))),data_ss2,'b.'); 
    hold on;  plot(.1*(0.5-rand(1,200))+zeros(1,200),[1+rand(1,200)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,8);
    threshold = 0.6; 
    data1 = [rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))+2*ones(1,numel(data_ss))),data_ss,'r.');
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))+2*ones(1,numel(data_ss2))),data_ss2,'b.'); 
    hold on;  plot(.1*(0.5-rand(1,200))+zeros(1,200),[1+rand(1,200)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,9);
    threshold = 1; 
    data1 = [rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))+2*ones(1,numel(data_ss))),data_ss,'r.'); 
    title('The "perfect" threshold');
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))+2*ones(1,numel(data_ss2))),data_ss2,'b.'); 
    hold on;  plot(.1*(0.5-rand(1,200)),[1+rand(1,200)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,10);
    threshold = 1.2; 
    data1 = 1+[rand(1,200)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,100)+2*ones(1,100)),[rand(1,100)],'b.'); 
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss)))+zeros(1,numel(data_ss)),data_ss,'r.'); 
    hold on;  plot(.1*(0.5-rand(1,numel(data_ss2)))+zeros(1,numel(data_ss2)),data_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,11);
    threshold = 1.6; 
    data1 = 1+[rand(1,200)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,100)+2*ones(1,100)),[rand(1,100)],'r.'); 
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss)))+zeros(1,numel(data_ss)),data_ss,'r.'); 
    hold on;  plot(.1*(0.5-rand(1,numel(data_ss2)))+zeros(1,numel(data_ss2)),data_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,12);
    threshold = 1.9; 
    data1 = 1+[rand(1,200)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    plot(.1*(0.5-rand(1,100)+2*ones(1,100)),[rand(1,100)],'r.'); 
    hold on;
    plot(.1*(0.5-rand(1,numel(data_ss)))+zeros(1,numel(data_ss)),data_ss,'r.'); 
    hold on;  plot(.1*(0.5-rand(1,numel(data_ss2)))+zeros(1,numel(data_ss2)),data_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    tbrate = [0.9,0.4,0,0,0,0];
    fbrate = [1,1,1,0.9,0.4,0.1];
    
    colors = hsv(6);
    
    subplot(3,6,13); 
    patch([0,1,1,0],[0,0,1,1],'r','edgecolor','none','facealpha',0); hold on; 
    for i = 1:6
        plot(fbrate(i),tbrate(i),'.','color',colors(i,:),'markersize',18);
    end
    text(0.3,0.3,sprintf('auROC~0'))
    xlabel('False blues'); ylabel('True blues'); box off;
    set(gca,'YTick',[0:.25:1]);
    ylim([0,1.5]);

    for i = 13:18
        subplot(3,6,i);
        ylim([0,1.5]);
        text(0.1,1.45,sprintf('True blues: %1.2f',tbrate(i-12)),'color',colors(i-12,:)); 
        text(0.1,1.3,sprintf('False blues: %1.2f',fbrate(i-12)),'color',colors(i-12,:)); 
        if i~=13; set(gca,'YTick',[],'XTick',[],'YColor','w','Xcolor','w'); end
    end 
    
    
    % AUC near 0.5
    figure('color','w'); 
    subplot(3,6,[1:3]); 
    plot([linspace(1,200,200)],[rand(1,200)],'r'); hold on;
    plot(300+linspace(1,100,100),[rand(1,100)],'k');
    plot([200+linspace(1,100,100)],[rand(1,100)],'b');
    title(sprintf('User can optionally specify the\nperiods of interest (red and blue)'))
    set(gca,'TickDir','out'); box off;
    subplot(3,6,[4:6]); plot(.1*(0.5-rand(1,200)+zeros(1,200)),[rand(1,100),rand(1,100)],'r.'); 
    hold on;  plot(.1*(0.5-rand(1,100))+.2*ones(1,100),[rand(1,100)],'b.'); xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;
    title(sprintf('A case with no separation\n between red and blue'))
    widebox = [0,0,1.2,.65]*1000;
    set(gcf,'position',widebox)

    subplot(3,6,7); 
    threshold = 0;
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    data2 = [rand(1,100)];
    data2_ss = data2(data2<threshold);
    data2_ss2 = data2(data2>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))),data_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))),data_ss2,'b.'); hold on;  
    plot(.1*(0.5-rand(1,numel(data2_ss))+2*ones(1,numel(data2_ss))),data2_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data2_ss2))+2*ones(1,numel(data2_ss2))),data2_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,8);
    threshold = 0.2; 
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    data2 = [rand(1,100)];
    data2_ss = data2(data2<threshold);
    data2_ss2 = data2(data2>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))),data_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))),data_ss2,'b.'); hold on;  
    plot(.1*(0.5-rand(1,numel(data2_ss))+2*ones(1,numel(data2_ss))),data2_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data2_ss2))+2*ones(1,numel(data2_ss2))),data2_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;
    
    subplot(3,6,9);
    threshold = 0.4; 
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    data2 = [rand(1,100)];
    data2_ss = data2(data2<threshold);
    data2_ss2 = data2(data2>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))),data_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))),data_ss2,'b.'); hold on;  
    plot(.1*(0.5-rand(1,numel(data2_ss))+2*ones(1,numel(data2_ss))),data2_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data2_ss2))+2*ones(1,numel(data2_ss2))),data2_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,10);
    threshold = 0.6;
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    data2 = [rand(1,100)];
    data2_ss = data2(data2<threshold);
    data2_ss2 = data2(data2>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))),data_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))),data_ss2,'b.'); hold on;  
    plot(.1*(0.5-rand(1,numel(data2_ss))+2*ones(1,numel(data2_ss))),data2_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data2_ss2))+2*ones(1,numel(data2_ss2))),data2_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;
    
    subplot(3,6,11);
    threshold = 0.8;
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    data2 = [rand(1,100)];
    data2_ss = data2(data2<threshold);
    data2_ss2 = data2(data2>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))),data_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))),data_ss2,'b.'); hold on;  
    plot(.1*(0.5-rand(1,numel(data2_ss))+2*ones(1,numel(data2_ss))),data2_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data2_ss2))+2*ones(1,numel(data2_ss2))),data2_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    subplot(3,6,12);
    threshold = 1;
    data1 = [rand(1,100),rand(1,100)];
    data_ss = data1(data1<threshold);
    data_ss2 = data1(data1>=threshold);
    data2 = [rand(1,100)];
    data2_ss = data2(data2<threshold);
    data2_ss2 = data2(data2>=threshold);
    plot(.1*(0.5-rand(1,numel(data_ss))),data_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data_ss2))),data_ss2,'b.'); hold on;  
    plot(.1*(0.5-rand(1,numel(data2_ss))+2*ones(1,numel(data2_ss))),data2_ss,'r.'); hold on;
    plot(.1*(0.5-rand(1,numel(data2_ss2))+2*ones(1,numel(data2_ss2))),data2_ss2,'b.'); 
    xlim([-.1,.3]);
    set(gca,'YTick',[],'YTickLabel',[],'YColor','w'); box off;

    tbrate = [1.0,0.8,0.6,0.4,0.2,0];
    fbrate = [1.0,0.8,0.6,0.4,0.2,0];

    colors = hsv(6);
    
    subplot(3,6,13); 
    patch([0,1,1,0],[0,0,1,0],'r','edgecolor','none','facealpha',0.2); hold on; 
    for i = 1:6
        plot(fbrate(i),tbrate(i),'.','color',colors(i,:),'markersize',18);
    end
    text(0.3,0.3,sprintf('auROC~0.5'))
    xlabel('False blues'); ylabel('True blues'); box off;
    set(gca,'YTick',[0:.25:1]);
    ylim([0,1.5]);

    for i = 13:18
        subplot(3,6,i);
        ylim([0,1.5]);
        text(0.1,1.45,sprintf('True blues: %1.2f',tbrate(i-12)),'color',colors(i-12,:)); 
        text(0.1,1.3,sprintf('False blues: %1.2f',fbrate(i-12)),'color',colors(i-12,:)); 
        if i~=13; set(gca,'YTick',[],'XTick',[],'YColor','w','Xcolor','w'); end
    end
    
       
end



end