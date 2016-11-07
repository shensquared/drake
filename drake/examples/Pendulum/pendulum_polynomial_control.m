function pendulum_polynomial_control()

doSolve=true;
doSimulation=true;
doPlotting=false;
deg_Y=0;
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

if doSolve
p = PendulumPlant();
n=3;
prog = spotsosprog;
x=msspoly('x',n);
% h=msspoly('h',n);
STATE=x;
A_x=[0, 0, (x(2)-1);0, 0, -x(1); -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
B=[0;0;1];
n_B=size(B,2);
x0=[0;0;0];
S0=x0;

prog = prog.withIndeterminate(STATE);
epsi=1e-4;

[prog,Q] = prog.newPSD(n);
% [prog,Y] = prog.newPSD(n);
% [prog,t] = prog.newPos(1);
% new_P=[X,t*eye(n);t*eye(n),Y];
V=x'*Q*x;
% V=STATE'*new_P*STATE;
prog = prog.withSOS(V-epsi*(STATE-S0)'*(STATE-S0));
% prog = prog.withSOS((x-x0)'*(x-x0)-V);
prog = prog.withEqs(subs(V,STATE,S0));  % V(0) = 0

m=monomials(x,0:deg_Y);
for i=1:n_B
	for j=1:n
        [prog,L_coeffs{i+j}] = prog.newFree(length(m));
        Y(i,j) = L_coeffs{i+j}'*m;
    end
end 

% % % % % % % % % % % 
Vdot=STATE'*([A_x*Q+Q*A_x'+B*Y+Y'*B'])*STATE;
deg_lambda = 4;
[prog,lambda] = prog.newFreePoly(monomials(x,0:deg_lambda));
prog = prog.withSOS(-Vdot-lambda*(x(1)^2+(x(2)-1)^2-1)-...
    -10*epsi*(STATE-S0)'*(STATE-S0));
prog = prog.withEqs(subs(Vdot,STATE,S0));
options = spot_sdp_default_options();
options.verbose = 0;
result = prog.minimize(-det(Q), solver,options);


Q=double(result.eval(Q))
for i=1:n_B
    for j=1:n
        L_coeffs{i+j} = double(result.eval(L_coeffs{i+j}));
        Y(i,j) = L_coeffs{i+j}'*m;
    end
end 

K=Y*inv(Q);
save('k.mat','K');
end

if doSimulation,
    pv = PendulumVisualizer();
    sys=MyPendulum2States;
    % random initial conditions
    initial_pos=randn(2,1);
    [ytraj,xtraj]=simulate(sys,[0,4],initial_pos);
    pv.playback(ytraj); 
    fnplt(xtraj,1);
end

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
end
