function [jrgeco,m] = SVDtoGECO(varargin)
filename = varargin{1};
load(filename)
sz = m.height/m.dsf;
if numel(varargin) == 2
    comps = varargin{2};
else
    comps = 1:500;
end

for c = 1:m.nLEDs
    Ctmp = C.(m.LEDs{c}); Stmp = S.(m.LEDs{c});
    tmp = zeros(sz^2,size(Stmp,2));
    tmp(m.nanidx,:) = Ctmp(:,comps)*Stmp(comps,:);
    tmp = tmp-min(tmp(:));
    data.(m.LEDs{c}) = reshape(uint16(tmp),[sz sz size(Stmp,2)]);
end
jrgeco = double(data.lime)./((double(data.red).^m.Dr).*(double(data.green).^m.Dg));
if isfield(m,'aux')
    m.rot = getrotf(m);
    m.bl = baselinefromrot(m.rot,100,2);
    m.bgGG = mean(jrgeco(:,:,m.bl),3);
else
    m.bl = 'no aux';
    m.bgGG = mean(jrgeco(:,:,30:130),3);
end
jrgeco = jrgeco./repmat(m.bgGG,[1 1 size(jrgeco,3)])-1;
jrgeco = reshape(jrgeco,[sz^2,size(jrgeco,3)]);
stdmap = std(jrgeco(:,1:100),[],2);
BW = stdmap~=0;
jrgeco(BW==0,:) = NaN;