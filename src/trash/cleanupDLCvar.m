function dataout = cleanupDLCvar(data)
data(find(data(:,3)<.99),1:2) = NaN;
data1 = hampel(sepblockfun(data(1:end-mod(size(data,1),3),1),[3 1],@mean));
data2 = hampel(sepblockfun(data(1:end-mod(size(data,1),3),2),[3 1],@mean));

dataout = [data1 data2];