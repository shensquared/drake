 function [K,L]=testingDOF()

% check solver dependency
checkDependency('spotless');
if checkDependency('mosek')
    solver = @spot_mosek;
else
    if ~checkDependency('sedumi')
        error('You will need to install Mosek or Sedumi to run this code');
    end
    solver = @spot_sedumi;
end
n=3;
prog = spotsosprog;
x=msspoly('x',n);
h=msspoly('h',n);
X=[x;h];
prog = prog.withIndeterminate(X);
x0=[0;-1;0];
h0=x0;
X0 = [x0;h0];

% % % % % % % candidate constant Lyapunov block element
% homogeneous Lyapunov function
deg_Q = 0;


% for i=1:n
%     for j=1:n
%         [prog,Q(i,j)]=prog.newFreePoly(monomials(h,deg_Q:deg_Q));
%     end
% end
% Q=Q+Q';



% or constant Q
[prog,Q] = prog.newPSD(n);
epsi=1e-9;

V=X'*[Q,zeros(n,n);zeros(n,n),Q]*X;
prog = prog.withSOS(V-epsi*(X-X0)'*(X-X0));
prog = prog.withEqs(subs(V,X,X0));  % V(0) = 0
% % % % % %
% % polynomial K, instrumented to be M, 
% only depends on the hat
deg_M = 3;
for i=1:n
    [prog,M(1,i)]=prog.newFreePoly(monomials(h,0:deg_M));
end

% m_monomials=prog.monomials(h,0:deg_M));
% [prog,M_coeffs]=prog.newFree(n*n*length(m_monomials)*length(m_monomials));
% M_coeffs=reshape(M_coeffs,n,n,length(m_monomials),length(m_monomials));
% [prog,M] = prog.newFreePoly(monomials(h,0:deg_M));


% % % % % % % % 
% % polynomial L, instrumented to be N, 
% % % % % % % % 
% depends on the measureable part of x and the hat


Y=[x(1);x(2);h(1);h(2);h(3)];

deg_N=3;
for i=1:n
    for j=1:n
        [prog, N(i,j)]=prog.newFreePoly(monomials(X,0:deg_N));
    end
end

p = PendulumPlant();
% V=X'*[Q,zeros(n,n);zeros(n,n),Q]*X;
% prog = prog.withEqs( subs(V,X,x0) );    % V(0) = 0

A_x=[0, x(3)/2, (x(2))/2;-x(3)/2, 0, -x(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
A_h=[0, h(3)/2, (h(2))/2;-h(3)/2, 0, -h(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];

% A_x=[0, x(3), 0;-x(3)/2, 0, -x(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
% A_h=[0, h(3), 0;-h(3)/2, 0, -h(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];

% A_x=[0, 0, (x(2));-x(3)/2, 0, -x(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
% A_h=[0, 0, (h(2));-h(3)/2, 0, -h(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];


% % % % % % % % % % % 

Vdot=X'*([A_x*Q,[0;0;1]*M;N,A_h*Q+[0;0;1]*M-N])*X;
deg_lambda = 3;
[prog,lambda] = prog.newFreePoly(monomials(x,0:deg_lambda));
[prog,lambda_hat] = prog.newFreePoly(monomials(h,0:deg_lambda));
prog = prog.withSOS(-Vdot-lambda*(x(1)^2+(x(2))^2-1)-lambda_hat*(h(1)^2+(h(2))^2-1)-epsi*(X-X0)'*(X-X0));

options = spot_sdp_default_options();
options.verbose = 1;
result = prog.minimize(-det(Q), solver,options);
Q=double((result.eval(Q)))
eig(Q)
K=result.eval(M)*inv(Q)
L=result.eval(N)*inv(Q)
save('DOFfeedback.mat','K','L')
% newstate=[sin(state(1));cos(state(1));state(2)];
% newstatehat=[sin(state(3));cos(state(3));state(4)];


% K=double(subs(result.eval(M),h,newstatehat))*inv(double(Q))*newstatehat
% K_x=double(subs(result.eval(M),h,newstate))*inv(double(Q))*newstate;


% L_gain_state=[double(subs(result.eval(N),X,[newstate;newstatehat]))]*inv(double(Q))*newstate;
% L_gain_state_hat=[double(subs(result.eval(N),X,[newstate;newstatehat]))]*inv(double(Q))*newstatehat;
% % translate L back to [q,qd] frame
% L_old=L_gain_state-L_gain_state_hat;

% L=[sqrt(L_old(1)^2+L_old(2)^2)*sign(state(1)-state(3));L_old(3)]
% L=[sqrt(L_gain_state(1)^2+L_gain_state(2)^2);L_gain_state(3)]-[sqrt(L_gain_state_hat(1)^2+L_gain_state_hat(2)^2);L_gain_state_hat(3)]

% L=[double(subs(result.eval(N),X,[newstate;newstatehat]))-double(subs(result.eval(N),X,[newstatehat]))]*inv(double(Q))*[newstate-newstatehat]
% L=double(subs(result.eval(N),X,state)/double(Q)*state;
end
