function inside = D(x) 
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: D_ex1_2.m
%--------------------------------------------------------------------------
% Description: Jump set
% Return 0 if outside of D, and 1 if inside D
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

    global T_t

    tauN = x(7);
    taud = x(8);
    q = x(9);
 
    tauT = x(18);
    tauM = x(19);
    tauS = x(20);
    qP = x(21);
    p = x(22);
    
    M_m = [x(23); x(24); x(25); x(26); x(27); x(28)];
    M_s = [x(29); x(30); x(31); x(32); x(33); x(34)];
     
%     
%     if (tauN <= 0 && q == 0) || (taud <= 0 && q == 1)
%         inside = 1;
%     else
%         inside = 0;
%     end
    
    if (tauT > T_t) || (tauN <= 0 && q == 0) || (taud <= 0 && q == 1)
        inside = 1;
    elseif (tauM < 0)
        if ((qP == 1) && (p == 1))
            inside = 1;
        end
        if ((qP == 1) && (p == 5))
            inside = 1;
        end
        if ((qP == 1) && (p == 4))
            inside = 1;
        end
    elseif (tauS < 0)
        if ((qP == 2) && (p == 3))
            inside = 1;
        end
        if ((qP == 2) && (p == 2))
            inside = 1;
        end    
    else
        inside = 0;
    end
    
end