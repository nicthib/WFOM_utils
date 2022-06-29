function TF = TF_within(mouse,day1,day2)
% This function outputs a TF for day 1 --> day 2 of the same mouse
m = makem;
m.dsf = 2; m.nrot = 2; m.BW = ones(256,256);
m.run = 'runB'; m.stim = '1';
m.outputs = 'g'; m.loadpct = [0 .01];
m.mouse = [mouse '_' day1];
fulldir = findmousefolder(m.mouse,m.run,m.stim);
[~,data] = LoadData_v2(fulldir,m);
tmp = max(data.green,[],3);
day1_im = uint8(round(255*tmp/max(tmp(:))));

m.mouse = [mouse '_' day2];
fulldir = findmousefolder(m.mouse,m.run,m.stim);
[~,data] = LoadData_v2(fulldir,m);
tmp = max(data.green,[],3);
day2_im = uint8(round(255*tmp/max(tmp(:))));

[optimizer, metric] = imregconfig('multimodal');
TF = imregister_nic(day1_im,day2_im,'similarity',optimizer,metric);
TF.im1 = day1_im;
TF.im2 = day2_im;


