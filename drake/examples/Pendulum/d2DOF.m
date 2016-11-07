function [K,L]=d2DOF()


% doPlotting=true;
doPlotting=false;
% doSimulation=false;
doSimulation=false;
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
A_x=[x(1), (x(1)); -x(1),0];
A_h=[h(1), (h(1)); -h(1),0];
B=[1;1];
n_B=size(B,2);
C=[1,0];
n_C=size(C,1);
x0=[0;0];
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
prog = prog.withSOS((1e4)*(x-x0)'*(x-x0)-half_V);

prog = prog.withEqs(subs(V,STATE,S0));  % V(0) = 0


% % % % % %
deg_K =2;
for i=1:n
    [prog,K_aux(n_B,i)]=prog.newFreePoly(monomials(h,0:deg_K));
end

% K_gain_var=K_aux(1,1)+K_aux(1,2);
% prog=prog.withSOS(1e3*h'*h-K_gain_var);
% prog=prog.withSOS(K_gain_var-1e3*h'*h);



deg_L1_aux=2;
observed=C*x;
for i=1:n
    for j=1:n_C
        [prog, L1_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L1_aux));
    end
end

deg_L2_aux=2;
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
result = prog.minimize(det(Y), solver,options);

X=double(result.eval(X))
eig(X)
Y=double(result.eval(Y))
eig(Y)
new_P=double(result.eval(new_P));
disp('K_aux');
K_aux=result.eval(K_aux)
disp('L1_aux')
L1_aux=result.eval(L1_aux)
disp('L2_aux')
L2_aux=result.eval(L2_aux)

disp('2nd SOS')
N=-Y;
eig(N)
M=(eye(n)-X*Y)*inv(N);
eig(M)

% new_prog = spotsosprog;
% new_prog = new_prog.withIndeterminate(STATE);
% [new_prog,M] = new_prog.newPSD(n);
% new_V=STATE'*([X,eye(n);eye(n),Y]-epsi*([X*X+M,X;X,eye(n)]))*STATE;
% new_prog=new_prog.withSOS(new_V);
% new_prog = new_prog.withEqs(subs(new_V,STATE,S0));
% new_Vdot=result.eval(Vdot);
% new_prog=new_prog.withSOS(-new_Vdot-epsi*STATE'*([X*X+M,X;X,eye(n)])*STATE);
% options.verbose = 0;
% sol=new_prog.minimize(-det(M),solver,options);
% MM=double(sol.eval(M));
% M=chol(MM);
% N=(inv(M)*(eye(n)-X*Y))';



pi_2=[eye(n),Y;zeros(n,n),N'];
eig(pi_2);
pi_1=[X,eye(n);M',zeros(n,n)];
eig(pi_1);
P=pi_2/pi_1;
eig(P);


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
K=result.eval(K_gain*h);
L1=result.eval(L1_gain*x);
L2=result.eval(L2_gain*h);

residule=L1-L2
subs(residule,[x1;h1;h3],[sin(2);sin(2);1])


% dmsubs(K,h,[sin(2);1]);

% K 1by3, L1  and L2 3by1
save('d2KL.mat','K','L1','L2');

total_V=STATE'*(P)*STATE;
real_vdot=STATE'*P*([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])*STATE;
pi_vdot=STATE'*(pi_1'*P*([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])*pi_1)*STATE;
difference=real_vdot-pi_vdot;
[Theta,ThetaDot] = meshgrid(-8:0.1:8, -8:0.25:8);  % for surfs

Vmesh = dmsubs(total_V,STATE,[Theta(:)';ThetaDot(:)';Theta(:)';ThetaDot(:)']);
Vmesh_min=min(Vmesh)
Vdotmesh=dmsubs(real_vdot,STATE,[Theta(:)';ThetaDot(:)';Theta(:)';ThetaDot(:)']);
vdot_max=max(Vdotmesh)
find(Vdotmesh>0)
pi_vdot_mesh=dmsubs(pi_vdot,STATE,[Theta(:)';ThetaDot(:)';Theta(:)';ThetaDot(:)']);
difference_mesh=dmsubs(difference,STATE,[Theta(:)';ThetaDot(:)';Theta(:)';ThetaDot(:)']);
if doPlotting,
    figure(1);
    % dmsubs(total_V,STATE,zeros(6,1))
    % dmsubs(real_vdot,STATE,zeros(6,1))
    % hhh=cos(Theta(:)')


    % Vmesh = dmsubs(total_V,STATE,[sin(Theta(:)');cos(Theta(:)');ThetaDot(:)';sin(Theta(:)');cos(Theta(:)');ThetaDot(:)']);
    % Vdotmesh=dmsubs(total_vdot,STATE,[sin(Theta(:)');cos(Theta(:)');ThetaDot(:)';sin(Theta(:)');cos(Theta(:)');ThetaDot(:)']);

    % Vmesh = dmsubs(total_V,STATE,[z(:)';y(:)';(:)';z(:)';y(:)';y3(:)']);
    subplot(1,2,1);
    surf(Theta,ThetaDot,reshape(Vmesh,size(Theta)));
    title('$$ V $$','interpreter','latex','fontsize',20) 
    xlabel('$$ theta $$','interpreter','latex','fontsize',15)
    ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

    subplot(1,2,2);
    surf(Theta,ThetaDot,reshape(Vdotmesh,size(Theta)));
    title('$$ real Vdot $$','interpreter','latex','fontsize',20) 
    xlabel('$$ theta $$','interpreter','latex','fontsize',15)
    ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

    figure(2);
    subplot(1,2,1);
    surf(Theta,ThetaDot,reshape(pi_vdot_mesh,size(Theta)));
    title('$$ pi_vdot $$','interpreter','latex','fontsize',20) 
    xlabel('$$ theta $$','interpreter','latex','fontsize',15)
    ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

    subplot(1,2,2);
    surf(Theta,ThetaDot,reshape(difference_mesh,size(Theta)));
    title('$$ difference $$','interpreter','latex','fontsize',20) 
    xlabel('$$ theta $$','interpreter','latex','fontsize',15)
    ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)


end



if doSimulation,
    sys=d2model;
    initial_pos=[1;1];
    initial_estimate=initial_pos;
    % initial_estimate=initial_pos+randn(2,1);
    initial_augmented=[initial_pos;initial_estimate];
    [ytraj,xtraj]=simulate(sys,[0,7],initial_augmented);
    fnplt(xtraj,1)
    fnplt(xtraj,2)
    fnplt(xtraj,3)
end



end
