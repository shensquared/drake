 function linearized_DOF()

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

p = PendulumPlant();
n=3;
prog = spotsosprog;
x=msspoly('x',n);
h=msspoly('h',n);
STATE=[x;h];
A_x=[0, 0, (-1);0, 0, 0; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
A_h=[0, 0, (-1);0, 0, 0; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
B=[0;0;1];
C=[1,0,0;0,1,0];
x0=[0;0;0];
h0=x0;
S0 = [x0;h0];

prog = prog.withIndeterminate(STATE);
epsi=1e-2;

[prog,X] = prog.newPSD(n);
[prog,Y] = prog.newPSD(n);
[prog,t] = prog.newPos(1);
new_P=[X,t*eye(n);t*eye(n),Y];
V=STATE'*new_P*STATE;
prog = prog.withSOS(V-epsi*(STATE-S0)'*(STATE-S0));
prog = prog.withEqs(subs(V,STATE,S0));  % V(0) = 0
% % % % % %

deg_K =0;
for i=1:n
    [prog,K_aux(1,i)]=prog.newFreePoly(monomials(h,0:deg_K));
end


deg_L1_aux=0;
observed=[x(1);x(2)];
observed_hat=[h(1);h(2)];
for i=1:n
    for j=1:2
        [prog, L1_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L1_aux));
    end
end

deg_L2_aux=0;
for i=1:n
    for j=1:3
        [prog, L2_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L2_aux));
    end
end


% % % % % % % % % % % 
Vdot=STATE'*([A_x*X+B*K_aux,A_x;L2_aux,Y*A_x+L1_aux*C])*STATE;
deg_lambda = 4;
[prog,lambda] = prog.newFreePoly(monomials(x,0:deg_lambda));
[prog,lambda_hat] = prog.newFreePoly(monomials(h,0:deg_lambda));
prog = prog.withSOS(-Vdot-lambda*(x(1)^2+(x(2)-1)^2-1)-lambda_hat*(h(1)^2+(h(2)-1)^2-1)-epsi*(STATE-S0)'*(STATE-S0));

options = spot_sdp_default_options();
options.verbose = 1;
result = prog.minimize(0, solver,options);


% P=double((result.eval(new_P)))
% det(P)
% eig(P)
X=double(result.eval(X));
eig(X)
Y=double(result.eval(Y));
eig(Y)

N=Y;
M=(X*Y-eye(n))*inv(Y);

K_aux=result.eval(K_aux);
K=K_aux*inv(M);



L1_aux=result.eval(L1_aux);
L1=inv(N)*L1_aux*C*(x);
% dmsubs(L1,[observed;h],[0;2;.1;2.1;5])

L2_aux=result.eval(L2_aux);




% [prog,N] = prog.newFree(n^2);
% N = reshape(N,n,n);
% A2=inv(N)*(L2_aux-(Y*A_x*X+L1_aux*C*(X-M)+Y*B*K_aux+N*B*K_aux))*inv(M)

% dmsubs(A_h(3,:)',[observed;h],[0;2;.1;2.1;5])
% dmsubs(A2(3,:)',[observed;h],[0;2;.1;2.1;5])

% save('partitionDOF_A.mat','K','L1','A2');

L2=inv(N)*((Y*A_x*X+L1_aux*C*X+Y*B*K_aux+N*A_h*M+N*B*K_aux)-L2_aux)*inv(M)*h;

residule=L1-L2

dmsubs(residule(1,:)',[observed;h],[0;2;0;2;1])
dmsubs(residule(2,:)',[observed;h],[0;2;0;2;1])
dmsubs(residule(3,:)',[observed;h],[0;2;0;2;1])

% dmsubs(residule(2,:)',[observed;h],[0;2;.1;2.1;3])
% dmsubs(residule(3,:)',[observed;h],[0;2;.1;2.1;3])

% K 1by3, L1 is 3by2, and L2 3by1
save('Linear_partitionDOF.mat','K','L1','L2');
end
