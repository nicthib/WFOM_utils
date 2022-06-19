function [data,m] = jrgeco_correction_2D(data,m,dbl)
F = data.lime./((data.red.^m.Dr).*(data.green.^m.Dg));
bl = dbl.lime./((dbl.red.^m.Dr).*(dbl.green.^m.Dg));
dF = F./bl; % F/F0
data.jrgeco = dF-1;