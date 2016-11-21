function DOF_nonunicircle()

doSimulation=false;
doPlotting=false;


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
A_x=[x(1), 0, (x(2)-1);x(1), 0, -x(1); -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
A_h=[h(1), 0, (h(2)-1);h(1), 0, -h(1); -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
B=[1;1;1];
C=[1,0,0;0,1,0];
x0=[0;0;0];
h0=x0;
S0 = [x0;h0];

prog = prog.withIndeterminate(STATE);
epsi=1e-4;

[prog,X] = prog.newPSD(n);
[prog,Y] = prog.newPSD(n);
[prog,t] = prog.newPos(1);
new_P=[X,t*eye(n);t*eye(n),Y];
half_V=x'*Y*x;
V=STATE'*new_P*STATE;
prog = prog.withSOS(V-epsi*(STATE-S0)'*(STATE-S0));
prog = prog.withSOS(-half_V+(1/epsi)*(x-x0)'*(x-x0));
prog = prog.withEqs(subs(V,STATE,S0));  % V(0) = 0
% % % % % %

deg_K =3;
for i=1:n
    [prog,K_aux(1,i)]=prog.newFreePoly(monomials(h,0:deg_K));
end


deg_L1_aux=3;
observed=[x(1);x(2)];
% observed_hat=[h(1);h(2)];
for i=1:n
    for j=1:2
        [prog, L1_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L1_aux));
    end
end

deg_L2_aux=3;
for i=1:n
    for j=1:3
        [prog, L2_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L2_aux));
    end
end


% % % % % % % % % % %

% half_Vdot=x'*A_x*X+B*K_aux
Vdot=STATE'*([A_x*X+B*K_aux,A_x;L2_aux,Y*A_x+L1_aux*C])*STATE;
% deg_lambda = 4;
% [prog,lambda] = prog.newFreePoly(monomials(x,0:deg_lambda));
% [prog,lambda_hat] = prog.newFreePoly(monomials(h,0:deg_lambda));
prog = prog.withSOS(-Vdot-epsi*(STATE-S0)'*(STATE-S0));

options = spot_sdp_default_options();
options.verbose = 1;
result = prog.minimize(det(Y)-det(X), solver,options);




X=double(result.eval(X))
eig(X);
Y=double(result.eval(Y))
eig(Y);
new_P=double(result.eval(new_P));
disp('K_aux');
K_aux=result.eval(K_aux)
disp('L1_aux')
L1_aux=result.eval(L1_aux)
disp('L2_aux')
L2_aux=result.eval(L2_aux)

disp('2nd SOS')
M=X;
eig(M)
N=(inv(M)*(eye(n)-X*Y))';
eig(N)


% new_prog = spotsosprog;
% new_prog = new_prog.withIndeterminate(STATE);
% [new_prog,M] = new_prog.newPSD(n);
% new_V=STATE'*([X,eye(3);eye(3),Y]-epsi*([X*X+M,X;X,eye(3)]))*STATE;
% new_prog=new_prog.withSOS(new_V);
% sol=new_prog.minimize(0,solver,options);
% MM=double(sol.eval(M));
% M=(chol(MM))'
% N=inv(M)*(eye(3)-X*Y)'

pi_2=[eye(n),Y;zeros(n,n),N']
eig_pi_2=eig(pi_2)
pi_1=[X,eye(n);M',zeros(n,n)]
eig_pi_1=eig(pi_1)
P=pi_2/pi_1;
eig_P=eig(P)


% only the gains 
disp('K_gain') 
K_gain=K_aux*inv(M')
disp('L1_gain')
L1_gain=inv(N)*L1_aux*C
disp('blah');
(Y*A_x*X+L1_aux*C*X+Y*B*K_aux+N*A_h*M'+N*B*K_aux)-L2_aux
disp('L2_gain');
L2_gain=inv(N)*((Y*A_x*X+L1_aux*C*X+Y*B*K_aux+N*A_h*M'+N*B*K_aux)-L2_aux)*inv(M')


% actual feedbacks
K=K_gain*h;
L1=L1_gain*x;
L2=L2_gain*h;




total_V=STATE'*(P)*STATE;
real_vdot=STATE'*P*([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])*STATE;
pi_vdot=STATE'*(pi_1'*P*([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])*pi_1)*STATE;
difference=real_vdot-pi_vdot;
[Theta,ThetaDot] = meshgrid(-pi:0.1:pi, -8:0.25:8);
    % hhh=cos(Theta(:)')

Vmesh = dmsubs(total_V,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
min(Vmesh)
Vdotmesh=dmsubs(real_vdot,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
max(Vdotmesh)
find(Vdotmesh>=0)



if doPlotting,
%%%%%%% plotting debug
total_V=STATE'*(P)*STATE;
total_vdot=STATE'*P*([A_x*X+B*K_aux,A_x;L2_aux,Y*A_x+L1_aux*C])*STATE;



[z,y] = meshgrid(-8:0.25:8,-8:0.25:8,-8:0.25:8);  % for surfs
figure(1);
% dmsubs(total_V,STATE,zeros(6,1))
% dmsubs(total_vdot,STATE,zeros(6,1))
Vmesh = dmsubs(total_V,STATE,[z(:)';y(:)';y3(:)';z(:)';y(:)';y3(:)']);
Vdotmesh=dmsubs(total_vdot,STATE,[z(:)';y(:)';y3(:)';z(:)';y(:)';y3(:)']);
subplot(1,2,1);
surf(z,y,reshape(Vmesh,size(z)));
title('$$ V $$','interpreter','latex','fontsize',20) 
xlabel('$$ x $$','interpreter','latex','fontsize',15)
ylabel('$$ y $$','interpreter','latex','fontsize',15)

subplot(1,2,2);
surf(z,y,reshape(Vdotmesh,size(z)));
title('$$ Vdot $$','interpreter','latex','fontsize',20) 
xlabel('$$ x $$','interpreter','latex','fontsize',15)
ylabel('$$ y $$','interpreter','latex','fontsize',15)

end
%%%%%%%% finishing plotting






% [prog,N] = prog.newFree(n^2);
% N = reshape(N,n,n);
% A2=inv(N)*(L2_aux-(Y*A_x*X+L1_aux*C*(X-M)+Y*B*K_aux+N*B*K_aux))*inv(M)

% dmsubs(A_h(3,:)',[observed;h],[0;2;.1;2.1;5])
% dmsubs(A2(3,:)',[observed;h],[0;2;.1;2.1;5])
% save('partitionDOF_A.mat','K','L1','A2');
% residule=L1-L2

% dmsubs(residule,[observed;h],zeros(5,1))
% dmsubs(residule,[observed;h],[sin(2);cos(2)+1;sin(2);cos(2)+1;1])

% K 1by3, L1  and L2 3by1
save('partitionDOF_nonuni.mat','K','L1','L2');


if doSimulation,
    sys=MyPendulum;
    initial_pos=[1;2;0];
    initial_estimate=initial_pos;
    % initial_estimate=initial_pos+randn(3,1);
    initial_augmented=[initial_pos;initial_estimate];
    [ytraj,xtraj]=simulate(sys,[0,7],initial_augmented);
end
end
