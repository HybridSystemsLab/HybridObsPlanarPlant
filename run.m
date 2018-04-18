%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: run_ex1_2.m
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system (bouncing ball)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00
%close all

global A M L T1 T2 T1_d T2_d q T_t d c

T1 = 0.2;
T2 = 1;
T1_d = 0;
T2_d = 0.2;
A=[0 1; -1 0];
M=[1 0];
L = [1.0097; 0.6015];
q = 0;
T_t = 3;
d = 0.5;
c = 0.5;

% plant initial conditions
z_0 = [10; 10;];

% plant estimate initial conditions
zhat_0 = [0; 0;];
zhat_nd_0 = [0; 0;];
zhat_p_0 = [0; 0;];

y_0 = M*z_0;
y_m0 = 0;

% tau initial condition

tauP0 = 0;
tauO0 = 10;
taud0 = 2*T2_d+1;
tauN0 = T1+(T2-T1).*rand(1,1);

tP_m0 = 0;
tP_mp0 = 0;

% Protocol Timer
tauT0 = 0;

% Master-Slave Timer
tauM0 = 0;

% Slave-Master Timer
tauS0 = 0;

% Slave Memory Buffer
M_m0 = zeros(6,1);

% Master Memory Buffer
M_s0 = zeros(6,1);

% Logic State Variables
q0 = 0;
qP0 = 0;
p0 = 1;

% 1588 state initial condition
w0 = [tauT0; tauM0; tauS0; qP0; p0; M_m0; M_s0;];

% state initial condition
x0 = [z_0; zhat_0; tauP0; tauO0; tauN0; taud0; q0; y_0; y_m0; tP_m0; zhat_nd_0; zhat_p_0; tP_mp0; w0];

% simulation horizon
TSPAN=[0 25];
JSPAN = [0 1000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.01);

% simulate
[t,j,x] = HyEQsolver( @f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

%%

% % plot solution
% figure(1) % state and estimate 1
% clf
% subplot(2,1,1), plotHarc(t,j,x(:,1));
% grid on
% ylabel('$z_1$','Interpreter','latex','FontSize',20)
% subplot(2,1,2), plotHarc(t,j,x(:,3));
% grid on
% ylabel('$\hat{z}_1$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
% 
% % plot solution
% figure(2) % state and estimate 2
% clf
% subplot(2,1,1), plotHarc(t,j,x(:,2));
% grid on
% ylabel('$z_2$','Interpreter','latex','FontSize',20)
% subplot(2,1,2), plotHarc(t,j,x(:,4));
% grid on
% ylabel('$\hat{z}_2$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
%%
% plot solution
figure(3) % error estimates (delay case)
clf
subplot(2,1,1), plotHarc(t,j,x(:,1) - x(:,3));
grid on
h = findobj(gca,'Type','line');
i = legend([h(1)],'$$\varepsilon_1$$');
set(i,'Interpreter','latex','FontSize',20)
%ylabel('$\varepsilon_1$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,2) - x(:,4));
grid on
h = findobj(gca,'Type','line');
i = legend([h(1)],'$$\varepsilon_2$$');
set(i,'Interpreter','latex','FontSize',20)
%ylabel('$\varepsilon_2$','Interpreter','latex','FontSize',20)
%xlabel('$t,j$','Interpreter','latex','FontSize',20)

%plot solution
figure(9) % error estimates (non-delay case)
clf
subplot(2,1,1), plotHarc(t,j,x(:,1) - x(:,13));
grid on
h = findobj(gca,'Type','line');
i = legend([h(1)],'$$\varepsilon_1$$');
set(i,'Interpreter','latex','FontSize',20)
%ylabel('$\varepsilon_1$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,2) - x(:,14));
grid on
h = findobj(gca,'Type','line');
i = legend([h(1)],'$$\varepsilon_2$$');
set(i,'Interpreter','latex','FontSize',20)
%ylabel('$\varepsilon_2$','Interpreter','latex','FontSize',20)
%xlabel('$t,j$','Interpreter','latex','FontSize',20)

%%

diff = abs(x(:,5) - x(:,6));

clear sync
for ii = 2:1:length(diff)
    if (abs(diff(ii)-diff(ii-1)) > 0.2)
        sync = t(ii);
    end
end
    
%%
modificatorF{1} = 'b';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1;
modificatorJ{1} = '--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1;
modificatorJ{4} = 'Marker';
modificatorJ{5} = '*';
modificatorJ{6} = 'MarkerEdgeColor';
modificatorJ{7} = 'b';
modificatorJ{8} = 'MarkerFaceColor';
modificatorJ{9} = 'b';
modificatorJ{10} = 'MarkerSize';
modificatorJ{11} = 5;

modificatorV{1} = 'r';
modificatorV{2} = 'LineWidth';
modificatorV{3} = 1;
modificatorM{1} = '--';
modificatorM{2} = 'LineWidth';
modificatorM{3} = 1;
modificatorM{4} = 'Marker';
modificatorM{5} = '*';
modificatorM{6} = 'MarkerEdgeColor';
modificatorM{7} = 'r';
modificatorM{8} = 'MarkerFaceColor';
modificatorM{9} = 'r';
modificatorM{10} = 'MarkerSize';
modificatorM{11} = 5;

%%
figure(3) % error estimates (delay case)
clf
subplot(2,1,1), plotHarc(t,j,x(:,1) - x(:,3),[],modificatorV,modificatorM);
hold on
subplot(2,1,1), plotHarc(t,j,x(:,1) - x(:,13),[],modificatorF,modificatorJ);
grid on
if exist('sync')
    vline(sync,'k','sync')
end
h = findobj(gca,'Type','line');
i = legend([h(130) h(181)],'$$\phi^{nom}$$','$$\phi^{\delta}$$');
%i = legend('$$\phi^{nom}$$','$$\phi^{\delta}$$');
%i = legend([h(35) h(2)],'$$norm_1(t,j)$$','$$norm_2(t,j)$$');
%jj = legend([h(3)],'$$\phi^{\delta}_{\varepsilon_1}$$');
set(i,'Interpreter','latex','FontSize',10)
%set(jj,'Interpreter','latex','FontSize',20)
ylabel('$\varepsilon_1$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,2) - x(:,4),[],modificatorV,modificatorM);
hold on
subplot(2,1,2), plotHarc(t,j,x(:,2) - x(:,14),[],modificatorF,modificatorJ);
grid on
if exist('sync')
    vline(sync,'k','sync')
end
% vline = refline([2292 0]);
% vline.Color = 'r';
%h = findobj(gca,'Type','line');
%i = legend('$$\phi^{nom}$$','$$\phi^{\delta}$$');
h = findobj(gca,'Type','line');
i = legend([h(130) h(181)],'$$\phi^{nom}$$','$$\phi^{\delta}$$');
%jj = legend([h(1)],'$$\phi^{\delta}_{\varepsilon_2}$$');
set(i,'Interpreter','latex','FontSize',10)
%set(jj,'Interpreter','latex','FontSize',20)
ylabel('$\varepsilon_2$','Interpreter','latex','FontSize',20)
xlabel('$t,j$','Interpreter','latex','FontSize',20)

%%


%%
% figure(5) % Plant and Observer timers
% clf
% plotHarc(t,j,x(:,5),[],modificatorV,modificatorM);
% hold on
% plotHarc(t,j,x(:,6),[],modificatorF,modificatorJ);
% grid on
% ylabel('$\tau_O$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)

% figure(4)
% clf
% plotHarc(t,j,diff,[],modificatorF,modificatorJ);
% grid on
% ylabel('$|\tau_P - \tau_O|$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
% 
% figure(5) % Network and Delay timers
% clf
% subplot(2,1,1), plotHarc(t,j,x(:,7));
% grid on
% ylabel('$\tau_N$','Interpreter','latex','FontSize',20)
% subplot(2,1,2), plotHarc(t,j,x(:,8));
% grid on
% ylabel('$\tau_{\delta}$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
% 
% modificatorF{1} = 'c';
% 
% figure(6) % Delay timer
% clf
% plotHarc(t,j,x(:,12));
% grid on
% ylabel('$\ell_{\tau_P}$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
% 
% 
% figure(7) % Y measurements 
% clf
% subplot(2,1,1), plotHarc(t,j,x(:,10));
% grid on
% ylabel('$y$','Interpreter','latex','FontSize',20)
% subplot(2,1,2), plotHarc(t,j,x(:,11));
% grid on
% ylabel('$y_{meas}$','Interpreter','latex','FontSize',20)
% xlabel('$t$','Interpreter','latex','FontSize',20)

%% comet plot

% plot solution
% figure(8) % error estimates (delay case)
% clf
% N = length(x(:,1));
% comet(x(N/2:N,1) - x(N/2:N,3),x(N/2:N,2) - x(N/2:N,4));
% %clf
% %plot(x(N/2:N,1) - x(N/2:N,3))

%%

e1 = x(:,1) - x(:,3);
e2 = x(:,2) - x(:,4);
e1_nd = x(:,1) - x(:,13);
e2_nd = x(:,2) - x(:,14);
e1_d_p = x(:,1) - x(:,15);
e2_d_p = x(:,1) - x(:,16);

e = [e1'; e2';];
e_nd = [e1_nd'; e2_nd';];
e_d_p = [e1_d_p'; e2_d_p';];

a = 0;
a_d = 0;
a_d2 = 0;
a_d_p = 0;
p = 0;
p_2 = 0;
p_const = 0;
p_p = 0;
norm_1 = 0;
norm_2 = 0;
%from the GainMatrixScript
P = [124.5712 -118.5413; -118.5413 226.1138];

for i = 1:length(e1_nd)
    a_d(i) = e(:,i)'*expm(A'*(x(i,7)))*P*expm(A*(x(i,7)))*e(:,i);
    a(i) = e_nd(:,i)'*expm(A'*(x(i,7)))*P*expm(A*(x(i,7)))*e_nd(:,i);
    a_d_p = e_d_p(:,i)'*expm(A'*(x(i,7)))*P*expm(A*(x(i,7)))*e_d_p(:,i);
    norm_1(i) = norm(e_nd(:,i));
    norm_2(i) = norm(e(:,i));
    if i > 1
        if j(i) - j(i-1) > 0
            if mod(j(i),2) == 1
                %p_const = (e_nd(:,i-1)' - e(:,i-1)')*expm(A'*(x(i-1,7)))*P*expm(A*(x(i-1,7)))*(e_nd(:,i-1) - e(:,i-1)) + (e(:,i)' - e(:,i)')*expm(A'*(x(i,7)))*P*expm(A*(x(i,7)))*(e(:,i) + e_nd(:,i));
                %p_const = (a(i-1)-a(i)) - (a_d(i-1) - a_d(i));
                Q1 = expm(A'*(x(i,7)))*P*expm(A*(x(i,7)));
                Q2 = expm(A'*(x(i-1,7)))*P*expm(A*(x(i-1,7)));
                %p_const = [e_nd(:,i)' e(:,i)']*[Q1 zeros(2); zeros(2) -Q1]*[e_nd(:,i); e(:,i);] + [e(:,i-1)' e_nd(:,i-1)']*[Q2 zeros(2); zeros(2) -Q2]*[e(:,i-1); e_nd(:,i-1);];
                %p_const = [e_nd(:,i)' e(:,i)']*[-Q1 zeros(2); zeros(2) Q1]*[e_nd(:,i); e(:,i);] + [e(:,i-1)' e_nd(:,i-1)']*[-Q2 zeros(2); zeros(2) Q2]*[e(:,i-1); e_nd(:,i-1);];
                %p_const2 = a_d(i-1) - a_d(i);
                p_const = a_d(i) - a(i);
            else
                %p(i) = 0;
                %p_const = a_d(i-1) - a_d(i);
                p_const = 0;
                %p_const2 = 0;
                %p_const = (a(i-1)-a(i)) - (a_d(i-1) - a_d(i));

            end
            %p_const = a_d(i-1) - a_d(i);
        end
    end
%     if mod(j(i),2) == 1
%         p(i) = p_const;
%     else
%         p(i) = 0;
%     end
    %p(i) = -1*p_const;
    p(i) = p_const;
    %p_2(i) = p_const2;
    a_d2(i) = a(i) + p(i);
end

view = [j'; x(:,7)'; p; e; a_d-a; e_nd; a_d; a;];

% modificatorF{1} = 'b';
% modificatorF{2} = 'LineWidth';
% modificatorF{3} = 3;
% modificatorJ{1} = '-.';
% modificatorJ{2} = 'LineWidth';
% modificatorJ{3} = 2;
% modificatorJ{4} = 'Marker';
% modificatorJ{5} = 'p';
% modificatorJ{6} = 'MarkerEdgeColor';
% modificatorJ{7} = 'r';
% modificatorJ{8} = 'MarkerFaceColor';
% modificatorJ{9} = 'b';
% modificatorJ{10} = 'MarkerSize';
% modificatorJ{11} = 6;
% 
% modificatorB{1} = 'r';
% modificatorB{2} = 'LineWidth';
% modificatorB{3} = 3;
% 
% %
% figure(10)
% subplot(2,1,1), plotHarc(t,j,a(:),[],modificatorF,modificatorJ);
% grid on
% ylabel('$\alpha(t,j)$','Interpreter','latex','FontSize',20)
% hold on
% subplot(2,1,2), plotHarc(t,j,a_d(:),[],modificatorF,modificatorJ);
% grid on
% % ylabel('$ \alpha_{\delta}(t,j)$','Interpreter','latex','FontSize',20)
% % subplot(4,1,3), plotHarc(t,j,x(:,7));
% % grid on
% % ylabel('$\tau_N$','Interpreter','latex','FontSize',20)
% % subplot(4,1,4), plotHarc(t,j,x(:,8));
% % grid on
% ylabel('$\alpha_{\delta}(t,j)$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)

%%

% figure(11)
% plotHarc(t,j,a_d(:),[],modificatorF,modificatorJ);
% hold on
% plotHarc(t,j,a(:),[],modificatorF,modificatorJ);
% grid on
% %ylabel('$$\alpha(t,j) , \alpha_{\delta}(t,j)$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',12)
% h = findobj(gca,'Type','line');
% i = legend([h(39) h(2)],'$$V \big (\phi_{\tilde{\mathcal{H}}}^{\delta}(t,j) \big )$$','$$V \big ( \phi_{\tilde{\mathcal{H}}}^N(t,j) \big )$$');
% set(i,'Interpreter','latex','FontSize',16)
%%


% figure(12)
% plotHarc(t,j,p(:),[],modificatorB,modificatorJ);
% hold on
% figure(6)
% plotHarc(t,j,abs(a_d(:) - a(:)),[],modificatorF,modificatorJ);
% grid on
% %ylabel('$$\alpha(t,j) , \alpha_{\delta}(t,j)$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
% h = findobj(gca,'Type','line');
% i = legend([h(35) h(2)],'$$\rho(t,j)$$','$$\alpha_{\delta}(t,j)$$');
% set(i,'Interpreter','latex','FontSize',20)
% 
% figure(13)
% plotHarc(t,j,a_d2(:),[],modificatorB,modificatorJ);
% hold on
% plotHarc(t,j,a(:),[],modificatorF,modificatorJ);
% grid on
% %ylabel('$$\alpha(t,j) , \alpha_{\delta}(t,j)$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
% h = findobj(gca,'Type','line');
% i = legend([h(35) h(2)],'$$\alpha_{\delta}(t,j)$$','$$\alpha(t,j)$$');
% set(i,'Interpreter','latex','FontSize',20)

%%
% figure(14)
% plotHarc(t,j,norm_1(:),[],modificatorB,modificatorJ);
% hold on
% plotHarc(t,j,norm_2(:),[],modificatorF,modificatorJ);
% grid on
% %ylabel('$$\alpha(t,j) , \alpha_{\delta}(t,j)$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
% h = findobj(gca,'Type','line');
% i = legend([h(161) h(2)],'$$||\varepsilon^{nom}||$$','$$||\varepsilon^{\delta}||$$');
% set(i,'Interpreter','latex','FontSize',16)
%%
% figure(9)
% plotHarc(t,j,e1(:),[],modificatorB,modificatorJ);
% hold on
% plotHarc(t,j,e1_nd(:),[],modificatorF,modificatorJ);
% grid on
% %ylabel('$$\alpha(t,j) , \alpha_{\delta}(t,j)$','Interpreter','latex','FontSize',20)
% xlabel('$t,j$','Interpreter','latex','FontSize',20)
% h = findobj(gca,'Type','line');
% i = legend([h(35) h(2)],'$$\varepsilon_1(t,j)$$','$$\varepsilon_1nd(t,j)$$');
% set(i,'Interpreter','latex','FontSize',20)

%% 

modificatorF{1} = 'b';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1;
modificatorJ{1} = '--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1;
modificatorJ{4} = 'Marker';
modificatorJ{5} = '*';
modificatorJ{6} = 'MarkerEdgeColor';
modificatorJ{7} = 'b';
modificatorJ{8} = 'MarkerFaceColor';
modificatorJ{9} = 'b';
modificatorJ{10} = 'MarkerSize';
modificatorJ{11} = 5;

figure(11)
plotHarc(t,j,norm_1(:),[],modificatorF,modificatorJ); %nominal
hold on
modificatorF{1} = 'r';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1;
modificatorJ{1} = '--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1;
modificatorJ{4} = 'Marker';
modificatorJ{5} = '*';
modificatorJ{6} = 'MarkerEdgeColor';
modificatorJ{7} = 'r';
modificatorJ{8} = 'MarkerFaceColor';
modificatorJ{9} = 'r';
modificatorJ{10} = 'MarkerSize';
modificatorJ{11} = 5;
plotHarc(t,j,norm_2(:),[],modificatorF,modificatorJ); %delay
grid on
if exist('sync')
    vline(sync,'k','sync')
end
ylabel('$||\varepsilon||$','Interpreter','latex','FontSize',20)
xlabel('$t,j$','Interpreter','latex','FontSize',20)
h = findobj(gca,'Type','line');
i = legend([h(161) h(2)],'$$\phi^{nom}$$','$$\phi^{\delta}$$');
set(i,'Interpreter','latex','FontSize',16)

%%
modificatorF{1} = 'b';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1;
modificatorJ{1} = '--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1;
modificatorJ{4} = 'Marker';
modificatorJ{5} = '*';
modificatorJ{6} = 'MarkerEdgeColor';
modificatorJ{7} = 'b';
modificatorJ{8} = 'MarkerFaceColor';
modificatorJ{9} = 'b';
modificatorJ{10} = 'MarkerSize';
modificatorJ{11} = 5;
figure(5) % Plant and Observer timers
clf
plotHarc(t,j,x(:,5),[],modificatorF,modificatorJ);
hold on
modificatorF{1} = 'r';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1;
modificatorJ{1} = '--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1;
modificatorJ{4} = 'Marker';
modificatorJ{5} = '*';
modificatorJ{6} = 'MarkerEdgeColor';
modificatorJ{7} = 'r';
modificatorJ{8} = 'MarkerFaceColor';
modificatorJ{9} = 'r';
modificatorJ{10} = 'MarkerSize';
modificatorJ{11} = 5;
plotHarc(t,j,x(:,6),[],modificatorF,modificatorJ);
grid on
%ylabel('$\tau_O$','Interpreter','latex','FontSize',20)
xlabel('$t,j$','Interpreter','latex','FontSize',20)
h = findobj(gca,'Type','line');
i = legend([h(161) h(2)],'$$\phi^{\delta}_{\tau_P}$$','$$\phi^{\delta}_{\tau_O}$$');
set(i,'Interpreter','latex','FontSize',16)