function d2linearDOF()



doPlotting=false;
doSimulation=true;
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

n=2;
prog = spotsosprog;
x=msspoly('x',n);
h=msspoly('h',n);
STATE=[x;h];
A_x=[3 5; 5,3];
A_h=[3 5; 5,3];
B=[1;1];
C=[1,0];
x0=[0;0];
h0=x0;
S0 = [x0;h0];

prog = prog.withIndeterminate(STATE);
epsi=1e-4;

[prog,X] = prog.newPSD(n);
[prog,Y] = prog.newPSD(n);
[prog,t] = prog.newPos(1);
new_P=[X,t*eye(n);t*eye(n),Y];

% half_V=x'*X*x;
V=STATE'*new_P*STATE;
prog = prog.withSOS(V-epsi*(STATE-S0)'*(STATE-S0));
% prog = prog.withSOS(half_V-epsi*(x-x0)'*(x-x0));
prog = prog.withEqs(subs(V,STATE,S0));  % V(0) = 0


% % % % % %
deg_K =1;
for i=1:n
    [prog,K_aux(1,i)]=prog.newFreePoly(monomials(h,0:deg_K));
end

deg_L1_aux=1;
observed=[x(1)];
for i=1:n
    for j=1:1
        [prog, L1_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L1_aux));
    end
end

deg_L2_aux=1;
for i=1:n
    for j=1:n
        [prog, L2_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L2_aux));
    end
end


% % % % % % % % % % % 
Vdot=STATE'*([A_x*X+B*K_aux,A_x;L2_aux,Y*A_x+L1_aux*C])*STATE;
% deg_lambda = 4;
% [prog,lambda] = prog.newFreePoly(monomials(x,0:deg_lambda));
% [prog,lambda_hat] = prog.newFreePoly(monomials(h,0:deg_lambda));
prog = prog.withSOS(-Vdot-epsi*(STATE-S0)'*(STATE-S0));

options = spot_sdp_default_options();
options.verbose = 0;
result = prog.minimize(det(X), solver,options);

X=double(result.eval(X));
eig(X);
Y=double(result.eval(Y));
eig(Y);
K_aux=result.eval(K_aux);
L1_aux=result.eval(L1_aux);
L2_aux=result.eval(L2_aux);

disp('2nd SOS')
% N=Y;
% eig(N)
% M=(eye(n)-X*Y)*inv(N);
% eig(M)

new_prog = spotsosprog;
new_prog = new_prog.withIndeterminate(STATE);
[new_prog,M] = new_prog.newPSD(n);

new_V=STATE'*([X,eye(n);eye(n),Y]-epsi*([X*X+M,X;X,eye(n)]))*STATE;
new_prog=new_prog.withSOS(new_V);

new_Vdot=result.eval(Vdot);
new_prog=new_prog.withSOS(-new_Vdot-epsi*STATE'*([X*X+M,X;X,eye(n)])*STATE);
options.verbose = 1;
sol=new_prog.minimize(0,solver,options);
MM=double(sol.eval(M));
M=chol(MM)
N=(inv(M)*(eye(n)-X*Y))'
pi_2=[eye(n),Y;zeros(n,n),N']
pi_1=[X,eye(n);M',zeros(n,n)]
P=pi_2/pi_1;
[eigvector,eigvalue]=eig(P)


% only the gains 
K_gain=K_aux*inv(M');
L1_gain=inv(N)*L1_aux*C;
L2_gain=inv(N)*((Y*A_x*X+L1_aux*C*X+Y*B*K_aux+N*A_h*M'+N*B*K_aux)-L2_aux)*inv(M');


% actual feedbacks

K=K_gain*h;
L1=L1_gain*x;
L2=L2_gain*h;

residule=L1-L2;
% dmsubs(residule,[observed;h],zeros(3,1))

dmsubs(residule,[observed;h],[sin(2);sin(2);1])
dmsubs(K,h,[sin(2);1])

% K 1by3, L1  and L2 3by1
save('d2KL_linear.mat','K','L1','L2');

if doPlotting,
total_V=STATE'*(P)*STATE;
total_vdot=STATE'*P*([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])*STATE;
[Theta,ThetaDot] = meshgrid(-8:0.25:8, -8:0.25:8);  % for surfs
figure(1);
% dmsubs(total_V,STATE,zeros(6,1))
% dmsubs(total_vdot,STATE,zeros(6,1))
% hhh=cos(Theta(:)');
% Vmesh = dmsubs(total_V,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
% Vdotmesh=dmsubs(total_vdot,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
Vmesh = dmsubs(total_V,STATE,[(Theta(:)');ThetaDot(:)';(Theta(:)');ThetaDot(:)']);
Vdotmesh=dmsubs(total_vdot,STATE,[(Theta(:)');ThetaDot(:)';(Theta(:)');ThetaDot(:)']);

% Vmesh = dmsubs(total_V,STATE,[z(:)';y(:)';(:)';z(:)';y(:)';y3(:)']);
subplot(1,2,1);
surf(Theta,ThetaDot,reshape(Vmesh,size(Theta)));

title('$$ V $$','interpreter','latex','fontsize',20) 
xlabel('$$ theta $$','interpreter','latex','fontsize',15)
ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

subplot(1,2,2);
surf(Theta,ThetaDot,reshape(Vdotmesh,size(Theta)));

title('$$ Vdot $$','interpreter','latex','fontsize',20) 
xlabel('$$ theta $$','interpreter','latex','fontsize',15)
ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)
end


if doSimulation,
    sys=d2model_linear;
    initial_pos=[1;2];
    initial_estimate=initial_pos;
    % initial_estimate=initial_pos+randn(2,1);
    initial_augmented=[initial_pos;initial_estimate];
    [ytraj,xtraj]=simulate(sys,[0,7],initial_augmented);
    fnplt(xtraj,1)
    fnplt(xtraj,2)
    fnplt(xtraj,3)
end



end
