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

% p = PendulumPlant();
n=2;
prog = spotsosprog;
x=msspoly('x',n);
h=msspoly('h',n);
STATE=[x;h];
A_x=[2,1;0,1]
B=[0;1];
C=[1,0];
n_B=size(B,2);
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
prog = prog.withSOS((1e5)*(x-x0)'*(x-x0)-half_V);
prog = prog.withEqs(subs(V,STATE,S0));  % V(0) = 0


% % % % % %
deg_C_hat =0;
for i=1:n
    [prog,C_hat(1,i)]=prog.newFreePoly(monomials(h,0:deg_C_hat));
end

deg_B_hat=0;
observed=C*x;
% observed_hat=[h(1);h(2)];
for i=1:n
    for j=1:n_C
        [prog, B_hat(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_B_hat));
    end
end

deg_A_hat=0;
for i=1:n
    for j=1:n
        [prog, A_hat(i,j)]=prog.newFreePoly(monomials([observed;h],0:deg_A_hat));
    end
end


% % % % % % % % % % % 
Vdot=STATE'*([A_x*X+B*C_hat,A_x;A_hat,Y*A_x+B_hat*C]+[A_x*X+B*C_hat,A_x;A_hat,Y*A_x+B_hat*C]')*STATE;
deg_lambda = 4;
% [prog,lambda] = prog.newFreePoly(monomials(x,0:deg_lambda));
% [prog,lambda_hat] = prog.newFreePoly(monomials(h,0:deg_lambda));
prog = prog.withSOS(-Vdot -10*epsi*(STATE-S0)'*(STATE-S0));
    % -lambda*(x(1)^2+(x(2)-1)^2-1)-lambda_hat*(h(1)^2+(h(2)-1)^2-1)...
   
% prog = prog.withEqs(subs(Vdot,STATE,S0));
options = spot_sdp_default_options();
options.verbose = 0;
result = prog.minimize(det(A_hat-A_x), solver,options);


X=double(result.eval(X))
eig(X);
Y=double(result.eval(Y))
eig(Y);
new_P=double(result.eval(new_P));
disp('C_hat');
C_hat=result.eval(C_hat)
disp('B_hat')
B_hat=result.eval(B_hat)
disp('A_hat')
A_hat=result.eval(A_hat)

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
eig_P=eig(P);


% only the gains 
disp('C_k') 
C_k=C_hat*inv(M')
disp('B_k')
B_k=inv(N)*B_hat
% disp('blah');
% (Y*A_x*X+B_hat*C*X+Y*B*C_hat+N*A_h*M'+N*B*C_hat)-A_hat
disp('A_k');
% A_k=inv(N)*((Y*A_x*X+B_hat*C*X+Y*B*C_hat+N*A_h*M'+N*B*C_hat)-A_hat)*inv(M')
A_k=inv(N)*(A_hat-N*B_hat*C*X-Y*B*C_hat)*inv(M')


% actual feedbacks
% K=K_gain*h;
% L1=L1_gain*x;
% L2=L2_gain*h;
% residule=L1-L2; 
% disp('residual');
% dmsubs(residule,[observed;h],[sin(0);cos(0)+1;sin(0);cos(0)+1;0])
% save('partitionDOF.mat','K','L1','L2');


% real_V=STATE'*(P)*STATE;
% real_vdot=STATE'*((P*[A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])+([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])'*P)*STATE;
% pi_times_realVdot=STATE'*(pi_1'*((P*[A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])+([A_x, B*K_gain;L1_gain,A_h+B*K_gain-L2_gain])'*P)*pi_1)*STATE;
% SOS_Vdot=result.eval(Vdot);
% [Theta,ThetaDot] = meshgrid(-pi:0.1:pi, -8:0.25:8);
% % hhh=cos(Theta(:)')
% Vmesh = dmsubs(real_V,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
% min(Vmesh)
% Vdotmesh=dmsubs(real_vdot,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
% max(Vdotmesh)
% pi_times_realVdot_mesh=dmsubs(pi_times_realVdot,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
% SOS_Vdot_mesh=dmsubs(SOS_Vdot,STATE,[sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)';sin(Theta(:)');cos(Theta(:)')+ones(1,4095);ThetaDot(:)']);
% ratio=pi_times_realVdot_mesh./SOS_Vdot_mesh;
% find(Vdotmesh>=0)


% if doPlotting,
%       % for surfs
%     figure(1);
%     subplot(1,2,1);
%     surf(Theta,ThetaDot,reshape(Vmesh,size(Theta)));
%     title('$$ real V $$','interpreter','latex','fontsize',20) 
%     xlabel('$$ theta $$','interpreter','latex','fontsize',15)
%     ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

%     subplot(1,2,2);
%     surf(Theta,ThetaDot,reshape(Vdotmesh,size(Theta)));
%     title('$$ real \dot{V} $$','interpreter','latex','fontsize',20) 
%     xlabel('$$ theta $$','interpreter','latex','fontsize',15)
%     ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

%     % figure(2);
%     % subplot(1,2,1);
%     % surf(Theta,ThetaDot,reshape(pi_times_realVdot_mesh,size(Theta)));
%     % title('$$ pi_times_realVdot $$','interpreter','latex','fontsize',20) 
%     % xlabel('$$ theta $$','interpreter','latex','fontsize',15)
%     % ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)

%     % subplot(1,2,2);
%     % surf(Theta,ThetaDot,reshape(difference_mesh,size(Theta)));
%     % title('$$ difference $$','interpreter','latex','fontsize',20) 
%     % xlabel('$$ theta $$','interpreter','latex','fontsize',15)
%     % ylabel('$$ thetadot $$','interpreter','latex','fontsize',15)


% end



% if doSimulation,
%     pv = PendulumVisualizer();
%     sys=MyPendulumCL;
%     initial_pos=[0;  12.9492];
%     initial_estimate=initial_pos;
% %     initial_estimate(2)=initial_pos(2)+randn(1,1);
%     initial_augmented=[initial_pos;initial_estimate]
%     [ytraj,xtraj]=simulate(sys,[0,50],[0;  12;0;12.9492]);
%     pv.playback(ytraj); 
%     fnplt(xtraj,1);
% end
end
