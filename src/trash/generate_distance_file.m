function generate_distance_file(filename,loadpath,savepath,win,cnt)
load(fullfile(loadpath,filename),'H')
cc = slcorr(H.jrgeco,1,win)';
D = pdist2(cc,cnt,'correlation')'; % label each full matrix using output centroids
save(fullfile(savepath,strrep(filename,'_H','_D')),'D','-v7.3')
