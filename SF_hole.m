function Truth = SF_hole(dir,Max)
%% This functions aims to randomly decide if we want to keep the input direction in the SF data set. The goal is to hole the Ewald sphere to gain realism.
%
% --- INPUTS --- 
% - 'dir' (1x3 vector) : the SF direction
% - 'Max' (float) : after this value, the probability to keep the direction
% is 0.
% ---------------
%
% --- OUTPUTS ---
% - 'Truth' (Boolean) : If True, the direction has been selected to be
% kept in data. If False, the contrary.
% ---------------
%
% YL.

    l = norm(dir,2);
    if l>Max
        Truth = false;
    else
        n=rand(1,3); % random draw for the 3 coordinates
        t = (1-exp(l/Max-1))/(1-exp(-1));
        if n(1)>t & n(2)>t & n(3)>t %In order for the hole to be made, all 3 coordinates need to be selected randomly at the same time
            Truth = false;
        else
            Truth = true;
        end
    end     
end