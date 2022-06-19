function generate_distance_file(filename,loadpath,savepath,opts)
load(fullfile(loadpath,filename),'H')
cc = slcorr(H.jrgeco,1,getwin(opts.ww,opts.sig));
d = pdist2(cc,opts.cnt,'correlation'); % label each full matrix using output centroids
save(fullfile(savepath,filename),'d','-v7.3')
