function Hnew = LPF(H,cutoffFreq)
ss = size(H);
zlen = 10000;
lpFilt = designfilt('lowpassiir','FilterOrder',20,...
'PassbandFrequency',cutoffFreq,'PassbandRipple',0.2,'SampleRate',19.9668);
hZ = [zeros(ss(1),zlen),H,zeros(ss(1),zlen)];
Hnew = [];
for n = 1:ss(1)
    Hnew(n,:) = filtfilt(lpFilt,hZ(n,:));
end
Hnew = Hnew(:,(zlen+1):(end-zlen));