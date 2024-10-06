%BUILD_QP Build a quadratic program for discrete linear optimal control
% Author: Daniel K. Mills
%
% Description: Build a quadratic program to solve the linear optimal
% control program:
% min e(:,L)'*QL*e(:,L) + sum_{k=1}^L e(:,k)'*Qk*e(:,k) + u(:,k)'*R*u(:,k)
% such that e(:,k) = x(:,k) - xr(:,k)
%           x(:,k+1) = A*x(:,k) + B*u(:,k)
%           x(:,1) = x0
%
% Inputs:
% L = Prediction horizon length
% QL = Terminal cost matrix (NxN)
% Q = Path state cost matrix (NxN) OR (NxNxL) for time-varying
% R = Path input cost matrix (MxM) OR (MxMxL) for time-varying
% A = Discrete linear state matrix (NxN) OR (NxNxL) for time-varying
% B = Discrete linear input matrix (NxM) OR (NxMxL) for time-varying
% x0 = Initial state (Nx1)
% xr = Reference state (Nx1) OR (Nx(L+1)) for time-varying leave empty '[]' if none
% xmin = Minimum state values (Nx1) OR (Nx(L+1)) for time-varying leave empty '[]' if none
% xmax = Maximum state values (Nx1) OR (Nx(L+1)) for time-varying leave empty '[]' if none
% umin = Minimum input values (Mx1) OR (NxL) for time-varying leave empty '[]' if none
% umax = Maximum input values (Mx1) OR (NxL) for time-varyingleave empty '[]' if none
%
% Outputs:
% H = Quadratic objective
% f = Linear objective
% Aeq = Linear constraint matrix Aeq*[x;u]=beq
% beq = Linear constraint constant Aeq*[x;u]=beq
% lb = Lower boundary lb <= [x;u]
% ub = Upper boundary [x; u] <= ub
function [H,f,Aeq,beq,lb,ub,M,N] = build_qp(L,QL,Q,R,A,B,x0,xr,xmin,xmax,umin,umax)
%% Prep
[N, M] = size(B(:,:,1));
if nargin < 8 || isempty(xr)
    xr = zeros(N,1);
end
if nargin < 9 || isempty(xmin)
    xmin = -Inf.*ones(N,1);
end
if nargin < 10 || isempty(xmax)
    xmax = Inf.*ones(N,1);
end
if nargin < 11 || isempty(umin)
    umin = -Inf.*ones(M,1);
end
if nargin < 12 || isempty(umax)
    umax = Inf.*ones(M,1);
end
%% Quadratic Program y = [x(0);...;x(L); u(0);...;u(L-1)]
if size(Q,3) == 1
    HQ = kron(speye(L),Q);
    fQtemp = repmat(-Q, L, 1);
else
    Qc = num2cell(Q,[1 2]);
    HQ = blkdiag(Qc{:});
    fQtemp = reshape(-Q,[],N);
end
if size(xr,2) == 1
    fQ = fQtemp*xr;
else
    fQ = zeros(N*L,1);
    ell = 1;
    for i = 1:N:(N*L)
        idx = i:(i+N-1);
        fQ(idx) = fQtemp(idx,:)*xr(:,ell);
        ell = ell + 1;
    end
end
if size(R,3) == 1
    HR = kron(speye(L),R);
else
    Rc = num2cell(R,[1 2]);
    HR = blkdiag(Rc{:});
end
H = blkdiag(HQ,QL,HR);  % Quadratic Objective
f = [fQ; -QL*xr(:,end); zeros(L*M, 1)];       % Linear Objective
Ax = -speye((L+1)*N);
if size(A,3) == 1
    Ax = Ax + kron(sparse(diag(ones(L,1),-1)),A);
else
    Ac = num2cell(A,[1 2]);
    Ax = Ax + [sparse(N,(L+1)*N); [sparse(blkdiag(Ac{:})) sparse(L*N,N)]];
end
if size(B,3) == 1
    Bu = kron([sparse(1,L);speye(L)],B);
else
    Bc = num2cell(B,[1 2]);
    Bu = [sparse(N,L*M); sparse(blkdiag(Bc{:}))];
end
%% Constraints
% Dynamics
Aeq = [Ax Bu];
beq = [-x0; zeros(L*N,1)];
% Box
if size(xmin,2) == 1
    xmin = repmat(xmin,L+1,1);
else
    xmin = xmin(:);
end
if size(xmax,2) == 1
    xmax = repmat(xmax,L+1,1);
else
    xmax = xmax(:);
end
if size(umin,2) == 1
    umin = repmat(umin,L,1);
else
    umin = umin(:);
end
if size(umax,2) == 1
    umax = repmat(umax,L,1);
else
    umax = umax(:);
end
lb = [xmin; umin];
ub = [xmax; umax];
end
