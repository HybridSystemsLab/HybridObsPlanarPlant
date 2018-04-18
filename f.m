function xdot = f(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: f_ex1_2.m
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system (bouncing ball)
% Description: Flow map
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00
    
    % state
    z = [x(1); x(2);];
    zhat = [x(3); x(4);];
    q = x(9);
    zhat_nd = [x(13); x(14);];
    zhat_dp = [x(15); x(16);];
    
    qP = x(21);
    
    M_m = [x(23); x(24); x(25); x(26); x(27); x(28)];
    M_s = [x(29); x(30); x(31); x(32); x(33); x(34)];
    
    global A M
    
    % differential equations
   
    if q == 1
        taudF = -1;
    else
        taudF = 0;
    end 
    
    % differential equations
    if qP == 1
        tauMdot = -1;
    else
        tauMdot = 0;
    end
    
    if qP == 2
        tauSdot = -1;
    else
        tauSdot = 0;
    end
    
    wdot = [1; tauMdot; tauSdot; 0; 0; zeros(6)*M_m; zeros(6)*M_s];
    
    xdot = [A*z; A*zhat; 1; 1; -1; taudF; 0; M*A*z; 0;0; A*zhat_nd; A*zhat_dp; 0; wdot];
    
end