function [n,f,Spectrum, Phase] = fastFourTrans(TimeSeries,Fs,ms)
% This functin calculates the fast fourier transfor of a one-dimensional time series.
% The signal sampling frequency should be fed to the function.

LoS = length(TimeSeries);       % Length of signal
n = 2^nextpow2(LoS);
f = Fs*(0:(n/2))/n;             % Range of trusted frequencies (anti-aliasing)

FFT = fft(TimeSeries,n);
Abs_FFT = abs(FFT/LoS);
angle_FFT = angle(FFT/LoS);
Spectrum = Abs_FFT(1:round(n/2)+1);       % Convert double-sided to single single-sided
Spectrum(2:end-1) = 2*Spectrum(2:end-1);    % Convert double-sided to single single-sided

Phase = angle_FFT(1:round(n/2)+1);       % Convert double-sided to single single-sided
Phase(2:end-1) = Phase(2:end-1);    % Convert double-sided to single single-sided
%figure;
%plot(f,Spectrum,'.')
%hold on
%plot(f,medfilt1(Spectrum,ms))
Spectrum = medfilt1(Spectrum,ms);
Phase = medfilt1(Phase,ms);
end