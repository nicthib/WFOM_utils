function wscore = getwhiskscore(whisk,win,thresh)

whisk = reshape(whisk,[numel(whisk) 1]);
whisk = padarray(whisk,[numel(win) 0],'pre');
for i = 1:numel(whisk)-numel(win)
    wscore(i) = sum(whisk(i:(i+numel(win)-1))>thresh)/numel(win);
    %wscore(i) = mean(whisk(i:(i+numel(win)-1)));
end
