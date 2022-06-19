function H = filter_H(H)

% 0.02-Hz high-pass filter to remove slow drifts (Ying Ma et al.,
% 2016:PNAS)
hpFilt = designfilt('highpassfir','StopbandFrequency',0.02, ...
    'PassbandFrequency',0.05,'PassbandRipple',0.5, ...
    'StopbandAttenuation',65,'DesignMethod','kaiserwin',...
    'SampleRate',19.9668);

% 2-Hz low-pass filter to reduce physiological noise (Ying et al.,
% 2016:PNAS)
lpFilt = designfilt('lowpassfir','PassbandFrequency',1.97, ...
    'StopbandFrequency',2,'PassbandRipple',0.5, ...
    'StopbandAttenuation',65,'DesignMethod','kaiserwin',...
    'SampleRate',19.9668);

% Pad
H.chbo = padarray(H.chbo,[0 100],0,'both');
H.chbr = padarray(H.chbr,[0 100],0,'both');
H.jrgeco = padarray(H.jrgeco,[0 100],0,'both');

for n = 1:size(H.jrgeco,1)
    H.chbo(n,:) = filtfilt(hpFilt,H.chbo(n,:));
    H.chbr(n,:) = filtfilt(hpFilt,H.chbr(n,:));
    H.chbo(n,:) = filtfilt(lpFilt,H.chbo(n,:));
    H.chbr(n,:) = filtfilt(lpFilt,H.chbr(n,:));
    H.jrgeco(n,:) = filtfilt(lpFilt,H.jrgeco(n,:));
end

H.chbo = H.chbo(:,101:end-100);
H.chbr = H.chbr(:,101:end-100);
H.jrgeco = H.jrgeco(:,101:end-100);