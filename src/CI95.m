% Returns 95% Confidence Interval for matrix y. Operates on columns.
function yCI95 = CI95(y)
y(isnan(y)) = 0;
N = size(y,1);
ySEM = std(y)/sqrt(N);
CI95 = tinv([0.975], N-1);
yCI95 = bsxfun(@times, ySEM, CI95(:));
