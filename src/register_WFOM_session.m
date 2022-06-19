
% TF = register_WFOM_session(FIXED,MOVING) performs an automated image
% registration between the images MOVING and FIXED, and outputs, TF, the
% affine transform needed to register MOVING to FIXED.
% Inputs;
% MOVING, FIXED: square WFOM images (preferrably from the green channel)
% Outputs:
% TF: affine2d variable that can be usied with imwarp() to register MOVING
% to FIXED.

function TF = register_WFOM_session(FIXED,MOVING)
sz = size(FIXED,1); % Get image size
cropidx = round(sz/5:4*sz/5); % Crop off first and last 20% of image
cropmask = nan(sz); cropmask(cropidx,cropidx) = 1; % Make cropmask

% Hist equalization
FIXED(~isnan(FIXED)) = histeq(uint16(FIXED(~isnan(FIXED))));
MOVING(~isnan(MOVING)) = histeq(uint16(MOVING(~isnan(MOVING))));

% Calculate imgradient and apply mask to MOVING and FIXED
FIXED_gradient = imgradient(FIXED.*cropmask);
% Normalize, convert to uint8 (improves imregister)
FIXED_gradient = uint8(255*FIXED_gradient/max(FIXED_gradient(:)));

% Repeat for MOVING
MOVING_gradient = imgradient(MOVING.*cropmask); 
MOVING_gradient = uint8(255*MOVING_gradient/max(MOVING_gradient(:)));

% Register MOVING to FIXED
[optimizer, metric] = imregconfig('monomodal');
TF = imregtform(MOVING_gradient,FIXED_gradient, 'similarity', optimizer, metric);

%imshowpair(FIXED,imwarp(MOVING,TF,'OutputView',imref2d(size(FIXED))))