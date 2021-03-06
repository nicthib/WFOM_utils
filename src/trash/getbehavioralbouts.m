function [running,onbouts,offbouts] = getbehavioralbouts(rot,thresh,pass1w,pass2w)
rot(rot < thresh) = 0;
% Pass 1: OR filter with width of pass1w frames
pass1 = zeros(size(rot));
for i = pass1w:numel(rot)
    if any(rot(i-pass1w+1:i) > 0)
        pass1(i) = 1;
    else
        pass1(i) = 0;
    end
end
% Larger pass of pass2w frames now that we have a clean signal
pass2 = zeros(size(pass1));
for i = pass2w:numel(pass2)
    if any(pass1(i-pass2w+1:i) > 0)
        pass2(i) = 1;
    else
        pass2(i) = 0;
    end
end
running = pass2;
if size(running,1) == 1
    running = running';
end
bouton = find(diff(running)==1);
boutoff = find(diff(running)==-1);
if numel(bouton) > numel(boutoff)
    boutoff = [boutoff; numel(rot)];
end
if numel(bouton) < numel(boutoff)
    bouton = [1; bouton];
end
for i = 1:numel(bouton)
    onbouts{i,1} = bouton(i):boutoff(i);
end
for i = 1:numel(bouton)-1
    offbouts{i,1} = boutoff(i):bouton(i+1);
end
