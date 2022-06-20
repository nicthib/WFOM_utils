function [chbo,chbr,chbt] = convert_Hb_2D(data1,data2,c1,c2,bl1,bl2,files)
%--------------------------------------------------------------------------
% 2D version of convert_Hb()
% convert takes in [blue,green] or [green,red] data sets and converts to
% HbO HbR and HbT 
%
% chb     : [HbR] output
% chbo    : [HbO] output (calculate HbT = HbO + HbR)
% data1   : first input data set (x,y,t)
% data2   : second input data set (x,y,t)
% c1      : color of data1 'b', 'g', or 'r'
% c2      : color of data2 'b','g', or 'r'
% filter  : 530 or 534
%
% Files to load - (files)
% dpffsMW_400to700nm.mat  : waves, dpff_488, dpff_530, dpff_630
% Hb_spectra.mat          : Hb, Hb02, lambda
% From https://omlc.org/spectra/hemoglobin/summary.html
% spectra_file            : spectra_green, spectra_blue, spectra_red
%--------------------------------------------------------------------------
fn = fieldnames(files);
for i = 1:numel(fn)
    load(files.(fn{i}));
end
spectra_blue(:,2) = spectra_blue(:,2)./max(spectra_blue(:,2));
spectra_green(:,2) = spectra_green(:,2)./max(spectra_green(:,2));
spectra_red(:,2) = spectra_red(:,2)./max(spectra_red(:,2));

x530 = spectra_green(:,1);
y530 = spectra_green(:,2);
x470 = spectra_blue(:,1);
y470 = spectra_blue(:,2);
x630 = spectra_red(:,1);
y630 = spectra_red(:,2);

Lolim = 400;
Hilim = 700;

splineyHb = spline([Lolim:2:Hilim],Hb(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
splineyHbO = spline([Lolim:2:Hilim],Hb02(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
spliney470 = spline(x470,y470,[Lolim:0.5:Hilim]);
spliney530 = spline(x530,y530,[Lolim:0.5:Hilim]);
spliney630 = spline(x630,y630,[Lolim:0.5:Hilim]);

splineydpff_488 = spline([Lolim:2:Hilim],dpff_488(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
splineydpff_530 = spline([Lolim:2:Hilim],dpff_530(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
splineydpff_630 = spline([Lolim:2:Hilim],dpff_630(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);

% ~ area under curves for E (do also for dpff)
% order is blue, green, red 
EHb(1) = sum((1/sum(spliney470))*spliney470.*splineyHb);
EHb(2) = sum((1/sum(spliney530))*spliney530.*splineyHb);
EHb(3) = sum((1/sum(spliney630))*spliney630.*splineyHb);
EHbO(1) = sum((1/sum(spliney470))*spliney470.*splineyHbO);
EHbO(2) = sum((1/sum(spliney530))*spliney530.*splineyHbO);
EHbO(3) = sum((1/sum(spliney630))*spliney630.*splineyHbO);

% still need to incorporate this into the DPF
DPF(1) = sum((1/sum(spliney470))*spliney470.*splineydpff_488);
DPF(2) = sum((1/sum(spliney530))*spliney530.*splineydpff_530);
DPF(3) = sum((1/sum(spliney630))*spliney630.*splineydpff_630);

% 1: blue, 2: green, 3: red
c1d = strfind('bgr',c1);
c2d = strfind('bgr',c2);

DPF1 = DPF(c1d);
DPF2 = DPF(c2d);
EHb1 = EHb(c1d);
EHb2 = EHb(c2d);
EHbO1 = EHbO(c1d);
EHbO2 = EHbO(c2d);

% make sure DPF and data are for same LED 
clear mua chbo chb chbt

%bl1 = mean(data1(:,:,round(baseinterval)),3);
%bl2 = mean(data2(:,:,round(baseinterval)),3);

mua(1,:,:)=-(1/DPF1)*log(data1./bl1);
mua(2,:,:)=-(1/DPF2)*log(data2./bl2);
clear data2 data1 DPF1 DPF2
chbo = squeeze((EHb2*mua(1,:,:)-EHb1*mua(2,:,:))/(EHbO1*EHb2-EHbO2*EHb1));
chbr = squeeze((EHbO2*mua(1,:,:)-EHbO1*mua(2,:,:))/(EHb1*EHbO2-EHb2*EHbO1));
chbt = chbo+chbr;
