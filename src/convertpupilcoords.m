function [out,perf] = convertpupilcoords(varargin)
data = varargin{1};
if numel(varargin) == 1
    thresh = 0;
else
    thresh = varargin{2};
end
out = [];
data.Eyelid = []; data.Mouse_Run = [];
data = rmfield(data,'Eyelid');
data = rmfield(data,'Mouse_Run');
parts_list = string(fieldnames(data));

% new approach
for i = 1:size(data.(parts_list(1)),1)
    x = []; y = [];
    for n = 1:(numel(parts_list))
        x(n) = data.(parts_list(n))(i,1);
        y(n) = data.(parts_list(n))(i,2);
        p(n) = data.(parts_list(n))(i,3);
    end
    x(p < thresh) = [];
    y(p < thresh) = [];
    if sum(p > thresh) == numel(parts_list)
        out(i,:) = CircleFitByPratt([x' y']);
    else
        out(i,1:3) = NaN;
    end
end
perf = sum(isnan(out(:,1)))/size(out,1);
