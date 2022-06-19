function [Q,bDr,bDg] = DrDg_brute_search(data,m,options)
arguments
    data struct
    m struct
    options.Dr (1,:) double = 0:.1:2 % Default Dr values to check
    options.Dg (1,:) double = 0:.05:.5 % Default Dg values to check
end

Q = zeros(numel(options.Dr),numel(options.Dg));
h = waitbar(0,'Calculating conversion scores...');
for i = 1:numel(options.Dr)
    for j = 1:numel(options.Dg)
        m.Dr = options.Dr(i);
        m.Dg = options.Dg(j);
        [~,m] = jrgeco_correction(data,m);
        Q(i,j) = m.jrgeco_conversion_score;
    end
    waitbar(i/numel(options.Dr),h)
end; close(h)

% Find minimum of 2D matrix. The is the approximate best score.
[~,idx] = min(Q(:));
[R,C] = ind2sub(size(Q),idx);
bDr = options.Dr(R);
bDg = options.Dg(C);
