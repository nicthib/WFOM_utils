function rotf = getrotf(m)
nFrames = round(m.movielength*m.framerate);
if ~isfield(m,'aux')
    m.aux = zeros(2,nFrames);
end
if iscell(m.aux) % converts m.aux cell to matrix
    tmpaux = [];
    disp('aux is bad!!')
    for d = 1:numel(m.aux)
        tmpaux(:,d) = m.aux{d}(1:2);
    end
    m.aux = double(tmpaux); clear tmpaux
end
rot = m.aux(2,:);
rot = rot(1:end-mod(numel(rot),3));
rot = sepblockfun(unwrap(rot,100),[1,3],'mean');
rot = [0 smooth(diff(rot))'];
rotf = rot;
