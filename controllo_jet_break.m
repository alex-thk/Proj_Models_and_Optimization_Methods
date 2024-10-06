%% Prep
clc;
clearvars;
close all;

% Simulation
t0 = 0;         % Initial time
tf = 1;         % Final time
L = 100;        % Horizon Length
Ts = 1/L;       % Discrete sample time
t = t0:Ts:tf;   % Discrete time vector

delta = 0.02; 
h = 1524; 
g = 9.81;
v0_x = 277.77; 
C = 0.45; 
% during the breaking phase alpha = 0  
A = 4; 
rho = 1.247; 
m = 10500; 
K = 10; 
d = 3500;
lambda = 2.7;  
beta = (0.5*C*A*rho)/m; 
gamma = K*delta*d;


A = [1 0 delta 0 0 0;
     0 1 0 0 0 0;
     -K*delta 0 (1-beta*delta) 0 1 0;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1
    ];


B= [delta;
    0;
    1; 
    0;
    0;
    0;
   ];

% Objective Function
QL = [1e9 0 0 0 0 0;   % Terminal cost matrix 
          0 1e-12 0 0 0 0;
          0 0 1e6 0 0 0;
          0 0 0 1e-12 0 0;
          0 0 0 0 1e-12 0;
          0 0 0 0 0 1e-12;
         ];       
    
Q = [1e3 0 0 0 0 0;     % Evolution cost matrix 
     0 1e-12 0 0 0 0;
     0 0 1e-3 0 0 0;
     0 0 0 1e-12 0 0; 
     0 0 0 0 1e-12 0;
     0 0 0 0 0 1e-12; 
    ];    

R = 1e10;           % Control cost matrix

% Box Constraints
xmax = [d+100; 0; 50; 0; gamma; 0];
xmin = [d; 0; 0; 0; gamma; 0];
umax = 10;
umin = -10;  

% Initial Conditions (of breaking phase)
x0 = [d; 0; 50; 0; gamma; 0];

% State reference (ideal trajectory) 
xr = [d+100; 0; 0; 0; gamma; 0];

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
subplot(3,1,1);
plot(t,x(1,:));
ylabel('Distance');
title('X Position of the jet on the runway');

subplot(3,1,2);
plot(t,x(3,:));
ylabel('Total X velocity');
 
subplot(3,1,3);
plot(t(1:L),u);
ylabel('Control X velocity');

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