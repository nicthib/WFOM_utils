% behav_states(rotf,ww,skipflag) assigns behavioral states defined by Bahar using rotf (rotary signal).
% Inputs:
% rotf: rotary signal (vector)
% ww: window width of epoch
% skipflag: boolean that returns either all state values, or only transition values.
% 
% Outputs:
% st: state vector (integer vector)
% idx: frame index for each values in st.
% rot_st: rotary values for each state (Nxww matrix, where N is number of states and ww is window width)

function [st,idx,rot_st] =  behav_states(rotf,ww,skipflag)
wr = (ww-1)/2;
rotf = rotf(:); % vectorize
st = []; idx = []; rot_st = [];
rotf = padarray(rotf,wr,'both');
rotf = smooth(rotf,20);
R = getrunningpulses(rotf);
RS_L = 40*20; % 40 seconds in
% N = not running, R = running
for i = RS_L+wr+1:numel(rotf)-ww*2
    epC = i-wr:i+wr;
    epP = epC - ww;
    epN = epC + ww;
    ep40 = epC - RS_L; % Particular to RS
    
    % State 1: Run onset. Defined as 50% N followed by, 50% R, with RS in
    % previous window and movement in next window
    if all(~R(epC(1:wr))) & all(R(epC(wr+2:end))) & all(R(epN)) & all(~R(epP))
        sttmp = 1; idxtmp = i;
    
    % State 2: Sustained R. Defined as 100% R
    elseif all(R(epC))
        sttmp = 2; idxtmp = i;
        
    % State 3: Run offset. Defined as the inverse of state 1.
    elseif  all(R(epC(1:wr))) & all(~R(epC(wr+2:end)))
        sttmp = 3; idxtmp = i;
        
    % State 4: Initial RS. Defined as no movement in current or next window, 
    % and movement in previous window
    elseif all(R(epP)) & all(~R(epC)) & all(~R(epN))
        sttmp = 4; idxtmp = i;
        
    % State 5: Sustained RS. Defined as no movement in previous 40 seconds
    % , current or next window.
    elseif all(~R(ep40)) & all(~R(epC)) & all(~R(epN))
        sttmp = 5; idxtmp = i;
    else
        sttmp = NaN; idxtmp = i;
    end

    if skipflag % skips repeated states, only looking for transitions
        if ~isempty(st) & st(end) == sttmp
            continue
        else
            if ~isnan(sttmp)
                st(end+1) = sttmp; 
                idx(end+1) = idxtmp;
                rot_st(end+1,:) = rotf(epC);
            end
        end
    else
        st(i) = sttmp;
    end
end


