function [on_frames,off_frames] = getrunonsetframes(runname,ploton)
load(runname,'rotf')
rz = getbehavioralbouts(rotf,1,20,10);
rz(1) = 0;
rz(end) = 0;
on_frames = find(diff(rz)==1);
off_frames = find(diff(rz)==-1);

if ploton
    figure
    plot(rotf)
    hold on
    plot(rz)
    scatter(on_frames,ones(numel(on_frames),1),'g','Filled')
    scatter(off_frames,ones(numel(off_frames),1),'r','Filled')
end