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

% plant estimate initial conditions for

% Observer with measurment delay
zhat_0 = [0; 0;];

% Observer with no measurement delay
zhat_nd_0 = [0; 0;];

% Auxiliary Observer variable for data analysis
zhat_p_0 = [0; 0;];

% Y measurements
y_0 = M*z_0;
y_m0 = 0;

% Initial condition for clocks and timers
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



%% Code to indicate when the clocks synchronize on the plots

diff = abs(x(:,5) - x(:,6));

clear sync
for ii = 2:1:length(diff)
    if (abs(diff(ii)-diff(ii-1)) > 0.2)
        sync = t(ii);
    end
end
    
%% Error estimate plots for both non-delay and delay Observers

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

figure(3)
clf
subplot(2,1,1), plotHarc(t,j,x(:,1) - x(:,3),[],modificatorV,modificatorM);
hold on
subplot(2,1,1), plotHarc(t,j,x(:,1) - x(:,13),[],modificatorF,modificatorJ);
grid on
if exist('sync')
    vline(sync,'k','sync')
end
h = findobj(gca,'Type','line');
i = legend([h(2) h(length(h)-1)],'$$\phi^{nom}$$','$$\phi^{\delta}$$');
set(i,'Interpreter','latex','FontSize',10)
ylabel('$\varepsilon_1$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,2) - x(:,4),[],modificatorV,modificatorM);
hold on
subplot(2,1,2), plotHarc(t,j,x(:,2) - x(:,14),[],modificatorF,modificatorJ);
grid on
if exist('sync')
    vline(sync,'k','sync')
end
h = findobj(gca,'Type','line');
i = legend([h(2) h(length(h)-1)],'$$\phi^{nom}$$','$$\phi^{\delta}$$');
set(i,'Interpreter','latex','FontSize',10)
ylabel('$\varepsilon_2$','Interpreter','latex','FontSize',20)
xlabel('$t,j$','Interpreter','latex','FontSize',20)

%% Code to generate the error norms

e1 = x(:,1) - x(:,3);
e2 = x(:,2) - x(:,4);
e1_nd = x(:,1) - x(:,13);
e2_nd = x(:,2) - x(:,14);
e1_d_p = x(:,1) - x(:,15);
e2_d_p = x(:,1) - x(:,16);

e = [e1'; e2';];
e_nd = [e1_nd'; e2_nd';];
e_d_p = [e1_d_p'; e2_d_p';];

norm_1 = 0;
norm_2 = 0;


for i = 1:length(e1_nd)
    norm_1(i) = norm(e_nd(:,i));
    norm_2(i) = norm(e(:,i));
end

%% Plot of the error norm for both delay and non-delay solutions

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
i = legend([h(length(h)-1) h(2)],'$$\phi^{nom}$$','$$\phi^{\delta}$$');
set(i,'Interpreter','latex','FontSize',16)

%% Plot of the Plant and Observer clocks

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
figure(5) 
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
xlabel('$t,j$','Interpreter','latex','FontSize',20)
h = findobj(gca,'Type','line');
i = legend([h(length(h)-1) h(2)],'$$\phi^{\delta}_{\tau_P}$$','$$\phi^{\delta}_{\tau_O}$$');
set(i,'Interpreter','latex','FontSize',16)