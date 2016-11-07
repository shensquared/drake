function [Q,K,L]=DOFpair()
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
epsi = 1e-9;
p = PendulumPlant();
x0 = [0;-1;0];
h0=x0;
X0=[x0;h0];
% % % % % % % candidate constant Lyapunov block element

[prog,Q] = prog.newPSD(n);
V=x'*Q*x;
prog.withSOS(V-epsi*(x-x0)'*(x-x0));
% % % % % %
% % polynomial K, instrumented to be M, 
% only depends on the hat
deg_M = 2;
for i=1:n
    [prog,M(1,i)]=prog.newFreePoly(monomials(h,0:deg_M));
end


% % % % % % % % 
% % polynomial L, instrumented to be N, 
% % % % % % % % 
% depends on both x and the hat
deg_N=2;
for i=1:n
    for j=1:n
        [prog, N(i,j)]=prog.newFreePoly(monomials(X,0:deg_N));
    end
end



A_x=[0, x(3)/2, x(2)/2;-x(3)/2, 0, -x(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
A_h=[0, h(3)/2, h(2)/2;-h(3)/2, 0, -h(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
B=[0;0;1];
% % % % % % % % % % % 
Vdot=X'*([A_x*Q,B*M;N,A_h*Q+B*M-N])*X;
% prog = prog.withEqs( subs(Q,X,x0) );  % V(0) = 0
deg_lambda = 2;
[prog,lambda] = prog.newFreePoly(monomials(X,0:deg_lambda));
[prog,lambda_hat] = prog.newFreePoly(monomials(X,0:deg_lambda));
prog = prog.withSOS(-Vdot-lambda*(x(1)^2+x(2)^2-1)-lambda_hat*(h(1)^2+h(2)^2-1));

options = spot_sdp_default_options();
options.verbose = 0;
result = prog.minimize(0, solver,options);
Q=double(result.eval(Q));
K=result.eval(M)*inv(Q)
L=result.eval(N)*inv(Q)
end