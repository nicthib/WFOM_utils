%% Initialize
clear
H_dir = '/local_mount/space/chronos/2/cmdata_analysis/RS_analysis/H_final'; % The base folder we are pulling data from
vars = {'H'}; % Variables required to include runs in our analysis
allrunnames = getallexps(H_dir,vars)'; % Gets all runs in the folder with the criteria above

%% Distance across subjects vs. window size
for i = 41:40:401
    %clearvars -EXCEPT allrunnames
    opts = [];
    opts.k = 12; opts.nreps = 100; % Some options
    opts.ww = i; opts.sig = 0; opts.dokmeans = 1;
    opts.skipfactor = 20;%(opts.ww-1)/2; % Some options
    runs = 'BCD'; mouse = 'cm12[5]_[3]';
    runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames)); % Parse out runs that end with B,C,D
    runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames)); % Parse out runs that are 'mouse'
    opts.variable = 'jrgeco';
    [~,c,~,dmin{1}] = getcorrstates(runnames,opts);
    
    mouse = 'cm12[6]_[3]';
    opts.dokmeans = 0; opts.cnt = c;
    runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,['run[' runs ']'])),allrunnames)); % Parse out runs that end with B,C,D
    runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames)); % Parse out runs that are 'mouse'
    [~,~,~,dmin{2}] = getcorrstates(runnames,opts);
    
    d1 = cell2mat(dmin{1}');
    d2 = cell2mat(dmin{2}');
    errorbar(i/20-.1,mean(d1(:)),std(d1(:)),'k')
    hold on
    errorbar(i/20+.1,mean(d2(:)),std(d2(:)),'r')
    drawnow
end
