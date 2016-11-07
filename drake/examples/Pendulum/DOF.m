function DOF()

doSimulation=true;
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
A_x=[0, 0, (x(2)-1);0, 0, -x(1); -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
A_h=[0, 0, (h(2)-1);0, 0, -h(1); -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
B=[0;0;1];
C=[1,0,0;0,1,0];
n_B=size(B,2);
n_C=size(C,1);
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
prog = prog.withSOS((1e5)*(x-x0)'*(x-x0)-half_V);
prog = prog.withEqs(subs(V,STATE,S0));  % V(0) = 0


% % % % % %
deg_K =2;
for i=1:n
    [prog,K_aux(1,i)]=prog.newFreePoly(monomials(h,0:deg_K));
end

deg_L1_aux=0;
observed=C*x;
% observed_hat=[h(1);h(2)];
for i=1:n
    for j=1:n_C
        [prog, L1_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L1_aux));
    end
end

deg_L2_aux=0;
for i=1:n
    for j=1:n
        [prog, L2_aux(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_L2_aux));
    end
end


% % % % % % % % % % % 
Vdot=STATE'*([A_x*X+B*K_aux,A_x;L2_aux,Y*A_x+L1_aux*C]+[A_x*X+B*K_aux,A_x;L2_aux,Y*A_x+L1_aux*C]')*STATE;
deg_lambda = 4;
[prog,lambda] = prog.newFreePoly(monomials(x,0:deg_lambda));
[prog,lambda_hat] = prog.newFreePoly(monomials(h,0:deg_lambda));
prog = prog.withSOS(-Vdot-lambda*(x(1)^2+(x(2)-1)^2-1)-lambda_hat*(h(1)^2+(h(2)-1)^2-1)...
    -10*epsi*(STATE-S0)'*(STATE-S0));
% prog = prog.withEqs(subs(Vdot,STATE,S0));
options = spot_sdp_default_options();
options.verbose = 0;
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


M=X;
% eig(M);
N=(inv(M)*(eye(n)-X*Y))';
% eig(N);


% disp('2nd SOS')
% new_prog = spotsosprog;
% new_prog = new_prog.withIndeterminate(STATE);
% [new_prog,M] = new_prog.newPSD(n);
% new_V=STATE'*([X,eye(n);eye(n),Y]-epsi*([X*X+M,X;X,eye(n)]))*STATE;
% new_prog=new_prog.withSOS(new_V);
% new_Vdot=result.eval(Vdot);
% new_prog=new_prog.withSOS(-new_Vdot-epsi*STATE'*([X*X+M,X;X,eye(n)])*STATE);
% options.verbose = 1;
% sol=new_prog.minimize(0,solver,options);
% MM=double(sol.eval(M));
% M=chol(MM)
% N=(inv(M)*(eye(3)-X*Y))'

pi_2=[eye(n),Y;zeros(n,n),N'];
eig_pi_2=eig(pi_2)
pi_1=[X,eye(n);M',zeros(n,n)];
eig_pi_1=eig(pi_1)
P=pi_2/pi_1;
eig_P=eig(P)


% only the gains 
disp('K_gain') 
K_gain=K_aux*inv(M')
disp('L1_gain')
L1_gain=inv(N)*L1_aux*C
% disp('blah');
% (Y*A_x*X+L1_aux*C*X+Y*B*K_aux+N*A_h*M'+N*B*K_aux)-L2_aux
disp('L2_gain');
L2_gain=inv(N)*((Y*A_x*X+L1_aux*C*X+Y*B*K_aux+N*A_h*M'+N*B*K_aux)-L2_aux)*inv(M')


% actual feedbacks
K=K_gain*h;
L1=L1_gain*x;
L2=L2_gain*h;
% residule=L1-L2; 
% disp('residual');
% dmsubs(residule,[observed;h],[sin(0);cos(0)+1;sin(0);cos(0)+1;0])
save('partitionDOF.mat','K','L1','L2');


real_V=STATE'*(P)*STATE;
real_vdot=STATE'*((P*[A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])+([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])'*P)*STATE;
pi_times_realVdot=STATE'*(pi_1'*((P*[A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])+([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])'*P)*pi_1)*STATE;
SOS_Vdot=result.eval(Vdot);
[Theta,ThetaDot] = meshgrid(-pi:0.1:pi, -8:0.25:8);
% hhh=cos(Theta(:)')
Vmesh = dmsubs(real_V,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
min(Vmesh)
Vdotmesh=dmsubs(real_vdot,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
max(Vdotmesh)
pi_times_realVdot_mesh=dmsubs(pi_times_realVdot,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
SOS_Vdot_mesh=dmsubs(SOS_Vdot,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
ratio=pi_times_realVdot_mesh./SOS_Vdot_mesh;
find(Vdotmesh>=0)


if doPlotting,
      % for surfs
    figure(1);
    subplot(1,2,1);
    surf(Theta,ThetaDot,reshape(Vmesh,size(Theta)));
    title('$$ real V $$','interpreter','latex','fontsize',20) 
    xlabel('$$ theta $$','interpreter','latex','fontsize',15)
    ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

    subplot(1,2,2);
    surf(Theta,ThetaDot,reshape(Vdotmesh,size(Theta)));
    title('$$ real \dot{V} $$','interpreter','latex','fontsize',20) 
    xlabel('$$ theta $$','interpreter','latex','fontsize',15)
    ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

    % figure(2);
    % subplot(1,2,1);
    % surf(Theta,ThetaDot,reshape(pi_times_realVdot_mesh,size(Theta)));
    % title('$$ pi_times_realVdot $$','interpreter','latex','fontsize',20) 
    % xlabel('$$ theta $$','interpreter','latex','fontsize',15)
    % ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

    % subplot(1,2,2);
    % surf(Theta,ThetaDot,reshape(difference_mesh,size(Theta)));
    % title('$$ difference $$','interpreter','latex','fontsize',20) 
    % xlabel('$$ theta $$','interpreter','latex','fontsize',15)
    % ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)


end



if doSimulation,
    pv = PendulumVisualizer();
    sys=MyPendulumCL;
    initial_pos=[0;  12.9492];
    initial_estimate=initial_pos;
%     initial_estimate(2)=initial_pos(2)+randn(1,1);
    initial_augmented=[initial_pos;initial_estimate]
    [ytraj,xtraj]=simulate(sys,[0,50],[0;  12;0;12.9492]);
    pv.playback(ytraj); 
    fnplt(xtraj,1);
end
end
