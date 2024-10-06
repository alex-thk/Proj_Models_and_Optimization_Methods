%% Prep
clc;
clearvars;
close all;

% Simulation
t0 = 0;         % Initial time
tf = 10;         % Final time
L = 1000;        % Horizon Length
Ts = 0.01;       % Discrete sample time
t = t0:Ts:tf;   % Discrete time vector

delta = 0.02; 
h = 1524; 
g = 9.81;
v0_x = 277.77; 
C = 0.45; 
% during the landing phase alpha is considered for simplicity = 3 
alpha = deg2rad(3);
Ah = 47.16; 
A = Ah*sin(alpha)+4;
rho = 1.247; 
m = 10500; 
K = 10; 
d = 3500;
lambda = 2.7;  
beta = (0.5*C*A*4*rho)/m; 
gamma = K*delta*d;


A = [1 0 delta 0 0 0;
     0 1 0 delta 0 0;
     0 0 (1-beta*delta) 0 1/L 0;
     0 0 0 1 delta 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1
    ];


B= [delta 0;
    0 delta;
    1 0; 
    0 1;
    0 0;
    0 0;
   ];

% Objective Function
QL = [1e9 0 0 0 0 0;   % Terminal cost matrix 
      0 1e9 0 0 0 0;
      0 0 1e12 0 0 0;
      0 0 0 1e9 0 0;
      0 0 0 0 1e-12 0;
      0 0 0 0 0 1e-12;
     ];    


Q = [1e9 0 0 0 0 0;     % Evolution cost matrix 
     0 1e-12 0 0 0 0;
     0 0 1e9 0 0 0;
     0 0 0 1e-12 0 0; 
     0 0 0 0 1e-12 0;
     0 0 0 0 0 1e-12; 
    ];    

R = [1 0;           % Control cost matrix 
     0 1
    ];           


% Box Constraints
xmax = [d; h; 350; 150; -10*lambda; -g];
xmin = [0; -h; -350; -150; -10*lambda; -g];
umax = [350; 150];
umin = [-350; -150];  

% Initial Conditions (of landing phase)
x0 = [0; h; 277.7; 0; -10*lambda; -g];

% State reference (ideal trajectory) 
xr = [d; 0; 50; 0; -10*lambda; -g];

%% Build the Quadratic Program
[H,f,Aeq,beq,lb,ub,M,N] = build_qp(L,QL,Q,R,A,B,x0,...
    xr,xmin,xmax,umin,umax);

%% MATLAB quadprog
y = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],...
    optimoptions('quadprog','Display','final-detailed',...
    'LinearSolver','sparse'));

% Results
[x,u] = parse_results(y,L,N,M);

figure;
subplot(2,1,1);
plot(t,x(1,:));
ylabel('Distance');
title('X Position of the jet while flying');

subplot(2,1,2);
plot(t,x(2,:));
ylabel('Altitude');
title('Y Position of the jet while flying');

figure; 
subplot(2,1,1);
plot(t,x(3,:));
ylabel('Total X velocity');

subplot(2,1,2);
plot(t,x(4,:));
ylabel('Total Y velocity'); 

figure; 
subplot(2,1,1);
plot(t(1:L),u(1,:));
ylabel('Control X velocity');

subplot(2,1,2);
plot(t(1:L),u(2,:));
ylabel('Control Y velocity');

% figure; 
% plot(t(1:L),u(1,:), LineWidth=1.5);
% ylabel('Control X velocity');
% title('Control X velocity in the first instant')
% xlim([0, 0.05])
% ylim([-5, 80])
% 
% figure; 
% plot(t,x(3,:), LineWidth=1.5);
% ylabel('Total X velocity');
% title('Total X velocity in the first instant')
% xlim([0, 0.05])
% ylim([250, 360])

% figure; 
% plot(t(1:L),u(1,:), LineWidth=1.5);
% ylabel('Control X velocity');
% title('Control X velocity in the last instant')
% xlim([9.95, 10])
% ylim([-5, 60])
% 
% figure; 
% plot(t(1:L),u(1,:), LineWidth=1.5);
% ylabel('Total X velocity');
% title('Total X velocity in the last instant')
% xlim([9.95, 10])
% ylim([0, 60])


%% Parse control states and inputs from quadratic program decision vector
function [x,u] = parse_results(y,L,N,M)
    x = zeros(N,L+1);
    temp = N*(L+1);
        for n = 1:N
            x(n,:) = y(n:N:temp);
        end
    u = zeros(M,L);
        for m = 1:M
            idx = (m:M:(M*L)) + temp;
            u(m,:) = y(idx);
        end
end