%% TF_a
% Transformation for cmw_x --> cmy_z
m = makem;
m.dsf = 2; m.nrot = 2; m.BW = ones(256,256);
m.run = ['runB']; m.stim = '1';
m.outputs = 'rgl'; m.loadpct = [0 .05]; 
m1 = 'cm125'; m2 = 'cm193';
% mt (mouse target) is always cm125
% ms is the mouse of interest
fulldir = findmousefolder([m1 '_1'],m.run,m.stim);
[~,mt] = LoadData_v2(fulldir,m);
fulldir = findmousefolder([m2 '_1'],m.run,m.stim);
[~,ms] = LoadData_v2(fulldir,m);

tmp = max(mt.green,[],3); mt.refim = uint8(round(255*tmp/max(tmp(:))));
tmp = max(ms.green,[],3); ms.refim = uint8(round(255*tmp/max(tmp(:))));
mt.p = cm125.p; ms.p = cm125.p;

%% Select registration points
% Here, we show the "fiducial" points on cm125's mask, and you are expected
% to adjust the othe fiducials on the right image to match. We use this to
% perform a cross mouse registration transformation matrix.
[mt.p,ms.p] = cpselect(mt.refim,ms.refim,mt.p,ms.p,'Wait',true);

%% Convert fiduciary points to masks
idx = convhull(mt.p(:,2),mt.p(:,1));
mt.regmask = uint8(poly2mask(mt.p(idx,1),mt.p(idx,2),256,256));
idx = convhull(ms.p(:,2),ms.p(:,1));
ms.regmask = uint8(poly2mask(ms.p(idx,1),ms.p(idx,2),256,256));

%% Register masks
[optimizer, metric] = imregconfig('monomodal');
varname = [m1 'to' m2];
TF_a.(varname) = imregister_nic(mt.regmask,ms.regmask,'affine',optimizer,metric);








