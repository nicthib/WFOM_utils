function rz = getrunningpulses(rotf,thresh,k,n)
% rz = getrunningpulses(rotf,thresh,k,n)
% This function generates pulses corresponding to running blocks.
% Inputs: rotf = rotary function; thresh = the threshold for generating
% pulses; k = the length of sliding window for movstd, moving standard
% deviation; n = the order of median filter; fig = on or off (if on plot
% the figure.

switch nargin
    case 0
    error('There is no input. Enter rotf')
    case 1
        thresh = 0.01;
        k = 20;
        n = 100;
    case 2
        k = 20;
        n = 100;
    case 3
        n = 100;
end


% calculate the derivative of rotf
rotf = rotf(:);
rotf_diff = diff(rotf);

% subtract the mode of rotf from rotf to make zero baseline; multiply the
% result by the derivative of rotf to make small changes smaller and large
% changes larger; this will allow to have smaller threshold
rotf_corrected = medfilt1((rotf - mode(rotf)),n) .* [rotf_diff;rotf_diff(end)];

% compute moving standard deviation 
rotf_corrected_std = movstd(rotf_corrected,k);

rotf_corrected_std(rotf_corrected_std > thresh) = 1;
rotf_corrected_std(rotf_corrected_std <= thresh) = 0;


% return rz
rz = rotf_corrected_std;

end