function [Gmag, Gdir] = imgradient_nic(varargin)

% Nic's version of imgradient:
% converts NaNs to 0, and I is forced to int16

narginchk(1,2);

[I, Gx, Gy, method] = parse_inputs(varargin{:});
% MY CODE HERE
I = int16(I);
I(isnan(I)) = 0;

% Compute directional gradients
if (isempty(I))     
    % Gx, Gy are given as inputs
    if ~isfloat(Gx)
        Gx = double(Gx);
    end
    if ~isfloat(Gy)
        Gy = double(Gy);
    end
    
else   
    % If Gx, Gy not given, compute them. For all others except Roberts
    % method, use IMGRADIENTXY to compute Gx and Gy. 
    if (strcmpi(method,'roberts'))        
        if ~isfloat(I)
            I = double(I);
        end
        Gx = imfilter(I,[1 0; 0 -1],'replicate');         
        Gy = imfilter(I,[0 1; -1 0],'replicate'); 
        
    else        
        [Gx, Gy] = imgradientxy(I,method);
        
    end
end

% Compute gradient magnitude
Gmag = hypot(Gx,Gy);

% Compute gradient direction
if (nargout > 1)
    if (strcmpi(method,'roberts'))
        
        Gdir = zeros(size(Gx));
        
        % For pixels with zero gradient (both Gx and Gy zero), Gdir is set
        % to 0. Compute direction only for pixels with non-zero gradient.
        xyNonZero = ~((Gx == 0) & (Gy == 0)); 
        Gdir(xyNonZero) = atan2(Gy(xyNonZero),-Gx(xyNonZero)) - (pi/4);
        Gdir(Gdir < -pi) = Gdir(Gdir < -pi) + 2*pi; % To account for the discontinuity at +-pi.
        
        Gdir = Gdir*180/pi; % Radians to degrees
        
    else
        
        Gdir = atan2(-Gy,Gx)*180/pi; % Radians to degrees
    end    
end

end
%======================================================================
function [I, Gx, Gy, method] = parse_inputs(varargin)

methodstrings = {'sobel','prewitt','roberts','centraldifference', ...
            'intermediatedifference'};
I = []; 
Gx = []; 
Gy = [];
method = 'sobel'; % Default method

if (nargin == 1)
    I = varargin{1};
    validateattributes(I,{'numeric','logical'},{'2d','nonsparse','real'}, ...
                       mfilename,'I',1);
        
else % (nargin == 2)
    if ischar(varargin{2}) || isstring(varargin{2})
        I = varargin{1};
        validateattributes(I,{'numeric','logical'},{'2d','nonsparse', ...
            'real'},mfilename,'I',1);
        validateattributes(varargin{2},{'char','string'}, ...
            {'scalartext'},mfilename,'METHOD',2);
        method = validatestring(varargin{2}, methodstrings, ...
            mfilename, 'METHOD', 2);
    else
        Gx = varargin{1};
        Gy = varargin{2}; 
        validateattributes(Gx,{'numeric','logical'},{'2d','nonsparse', ...
                           'real'},mfilename,'Gx',1);
        validateattributes(Gy,{'numeric','logical'},{'2d','nonsparse', ...
                           'real'},mfilename,'Gy',2);
        if (~isequal(size(Gx),size(Gy)))
            error(message('images:validate:unequalSizeMatrices','Gx','Gy'));
        end
    end
         
end

end
%----------------------------------------------------------------------