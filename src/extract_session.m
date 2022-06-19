function epochs = extract_session(rotf)
R = getrunningpulses(rotf);
R(1) = 0; R(end) = 0;
ons = []; ofs = [];
for i=1:numel(R)-1
    if R(i) == 1 && R(i+1) == 0
        ofs(end+1) = i;
    end
    if R(i) == 0 && R(i+1) == 1
        ons(end+1) = i;
    end
end

epochs = [];
if isempty(ons)
    return
end
epochs(1,:) = [1 ons(1)];
for i=1:numel(ons)-1
    epochs(end+1,:) = [ofs(i) ons(i+1)];
end
epochs(end+1,:) = [ofs(end) numel(R)];
