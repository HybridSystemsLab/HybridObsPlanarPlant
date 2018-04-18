function xplus = g(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: g_ex1_2.m
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system (bouncing ball)
% Description: Jump map
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

    % state

    global A M L T1 T2 T1_d T2_d c d T_t
    
    z = [x(1); x(2);];
    zhat = [x(3); x(4);];
    tauP = x(5);
    tauO = x(6);
    tauN = x(7);
    taud = x(8);
    q    = x(9);
    y    = x(10);
    y_m  = x(11);
    tP_m = x(12);
    zhat_nd = [x(13); x(14);];
    zhat_p = [x(15); x(16);];
    tP_mp = x(17);

    tauT = x(18);
    tauM = x(19);
    tauS = x(20);
    qP = x(21);
    p = x(22);
    
    M_m = [x(23); x(24); x(25); x(26); x(27); x(28)];
    M_s = [x(29); x(30); x(31); x(32); x(33); x(34)];
     
%     % 1588 state initial condition
% w0 = [tauT0; tauM0; tauS0; qP0; M_m0; M_s0; p0];
% 
% % state initial condition
% x0 = [z_0; zhat_0; tauP0; tauO0; tauN0; taud0; q0; y_0; y_m0; tP_m0; zhat_nd_0; zhat_p_0; tP_mp0; w0];

%     C1 = 10;
%     C2 = 1;
    
    G = (1 - q) + 2*q;
    
%     tauT_plus = tauT;
%     tauM_plus = tauM;
%     tauS_plus = tauS;
%     
%     qP_plus = qP;
%     p_plus = p;
% 
%     M_m_plus = M_m;
%     M_s_plus = M_s;
    
    if (tauN <= 0) && (G == 1)
        zhat_plus = zhat;
        zhat_p_plus = zhat_p;
        taud_plus = T1_d+(T2_d-T1_d).*rand(1,1);
        tauN_plus = T1+(T2-T1).*rand(1,1);
        y_m       = y;
        tP_m      = tauP;
        tP_mp     = tauP; %+ C1*exp(C2*tauP);
        q         = 1;
        zhat_nd_plus   = zhat_nd + L*(y_m - M*zhat_nd);
        
        %protocol
        tauT_plus = tauT;
        tauM_plus = tauM;
        tauS_plus = tauS;
        
        qP_plus = qP;
        p_plus = p;
        
        M_m_plus = M_m;
        M_s_plus = M_s;
        
        tauP_plus = tauP;
        tauO_plus = tauO;
        
    elseif (taud <= 0) && (G == 2)
        
        %non-perturbed
        del = tauO - tP_m;
        zhat_del  = expm(-A*del)*zhat;
        zhat_rev  = zhat_del + L*(y_m - M*zhat_del);
        zhat_plus = expm(A*del)*zhat_rev;
        
        %perturbed
        del_p = tauO - tP_mp;
        zhat_p_del  = expm(-A*del_p)*zhat_p;
        zhat_p_rev  = zhat_p_del + L*(y_m - M*zhat_p_del);
        zhat_p_plus = expm(A*del_p)*zhat_p_rev;
        
        taud_plus = 2*T2_d + 1;
        tauN_plus = tauN;
        q = 0;
        zhat_nd_plus = zhat_nd;
        
        %protocol
        tauT_plus = tauT;
        tauM_plus = tauM;
        tauS_plus = tauS;
        
        qP_plus = qP;
        p_plus = p;
        
        M_m_plus = M_m;
        M_s_plus = M_s;
        
        tauP_plus = tauP;
        tauO_plus = tauO; 
    
    elseif tauT >= T_t
        
        tauT_plus = 0;
        tauM_plus = d;
        tauS_plus = tauS;
        
        qP_plus = 1;
        p_plus = 1;
        
        M_m_plus = [tauP; 0; 0; 0; 0; 0];
        M_s_plus = M_s;
        
        tauP_plus = tauP;
        tauO_plus = tauO;        
    
        zhat_plus = zhat; 
        tauP_plus = tauP; 
        tauO_plus = tauO; 
        tauN_plus = tauN; 
        taud_plus = taud; 

        zhat_nd_plus = zhat_nd; 
        zhat_p_plus = zhat_p; 
    
    elseif tauM < 0
        
        zhat_plus = zhat; 
        tauP_plus = tauP; 
        tauO_plus = tauO; 
        tauN_plus = tauN; 
        taud_plus = taud; 
        zhat_nd_plus = zhat_nd; 
        zhat_p_plus = zhat_p;
        
        if qP == 1
            %Step 5
            if p == 4
                tauT_plus = tauT;
                tauM_plus = d;
                tauS_plus = 0;

                qP_plus = 1;
                p_plus = p + 1;

                M_m_plus = [M_m(1); M_m(2); M_m(3); M_m(4); tauP; 0];
                M_s_plus = M_s;

                tauP_plus = tauP;
                tauO_plus = tauO;
            %Step 2
            elseif p == 1
                tauT_plus = tauT;
                tauM_plus = 0;
                tauS_plus = c;

                qP_plus = 2;
                p_plus = p + 1;

                M_m_plus = M_m;
                M_s_plus = [M_m(1); tauO; 0; 0; 0; 0];

                tauP_plus = tauP;
                tauO_plus = tauO;
            %step 6
            elseif p == 5
                tauT_plus = tauT;
                tauM_plus = 0;
                tauS_plus = 0;

                qP_plus = 0;
                p_plus = 0;

                M_m_plus = M_m;
                M_s_plus = [M_m(1); M_m(2); M_m(3); M_m(4); M_m(5); tauO];

                tauP_plus = tauP;
                
                
                if (tauO > tauP + 0.01) || (tauO < tauP - 0.01)                  
                    tauO_plus = tauO - (M_m(2)-M_m(1)-M_m(4)+M_m(3))/2;
                else
                    tauO_plus = tauO;
                end
            end
        end
    
    
    
    elseif tauS < 0
        
        zhat_plus = zhat; 
        tauP_plus = tauP; 
        tauO_plus = tauO; 
        tauN_plus = tauN; 
        taud_plus = taud; 
        zhat_nd_plus = zhat_nd; 
        zhat_p_plus = zhat_p;
        
        if qP == 2
            %Step 3
            if p == 2
                tauT_plus = tauT;
                tauM_plus = 0;
                tauS_plus = d;

                qP_plus = 2;
                p_plus = p + 1;

                M_m_plus = M_m;
                M_s_plus = [M_s(1); M_s(2); tauO; 0; 0; 0];

                tauP_plus = tauP;
                tauO_plus = tauO;
            % Step 4
            elseif p == 3
                tauT_plus = tauT;
                tauM_plus = c;
                tauS_plus = 0;

                qP_plus = 1;
                p_plus = p + 1;

                M_m_plus = [M_s(1); M_s(2); M_s(3); tauP; 0; 0];
                M_s_plus = M_s;

                tauP_plus = tauP;
                tauO_plus = tauO;
            end
        end
    end
    
    wplus = [tauT_plus; tauM_plus; tauS_plus; qP_plus; p_plus; M_m_plus; M_s_plus;];
    xplus = [z; zhat_plus; tauP_plus; tauO_plus; tauN_plus; taud_plus; q; y; y_m; tP_m; zhat_nd_plus; zhat_p_plus; tP_mp; wplus];
    
end