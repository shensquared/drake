function plotingDOF()

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
x0 = [0;-1;0;0;-1;0];


prog = prog.withIndeterminate(X);

% homogeneous Lyapunov function
% deg_Q = 4;
% for i=1:n
%     for j=1:n
%         [prog,Q(i,j)]=prog.newFreePoly(monomials(X,deg_Q:deg_Q));
%     end
% end
% Q=Q+Q';
[prog,Q] = prog.newPSD(n);
epsi=1e-8;
V=X'*[Q,zeros(n,n);zeros(n,n),Q]*X;
prog = prog.withSOS(V-epsi*(X-x0)'*(X-x0));
% prog = prog.withEqs(subs(V,X,x0));  % V(0) = 0

deg_M = 2;
for i=1:n
    [prog,M(1,i)]=prog.newFreePoly(monomials(h,0:deg_M));
end

deg_N=2;
for i=1:n
    for j=1:n
        [prog, N(i,j)]=prog.newFreePoly(monomials(X,0:deg_N));
    end
end

p = PendulumPlant();


% prog = prog.withEqs( subs(V,X,x0) );    % V(0) = 0
A_x=[0, x(3)/2, (x(2))/2;-x(3)/2, 0, -x(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];
A_h=[0, h(3)/2, (h(2))/2;-h(3)/2, 0, -h(1)/2; -p.g/p.l, 0, -p.b/(p.m*p.l.^2)];

% % % % % % % % % % % 

Vdot=X'*([A_x*Q,[0;0;1]*M;N,A_h*Q+[0;0;1]*M-N])*X;
deg_lambda = 2;
[prog,lambda] = prog.newFreePoly(monomials(X,0:deg_lambda));
[prog,lambda_hat] = prog.newFreePoly(monomials(X,0:deg_lambda));
prog = prog.withSOS(-Vdot-lambda*(x(1)^2+(x(2))^2-1)-lambda_hat*(h(1)^2+(h(2))^2-1)-epsi*(X-x0)'*(X-x0)*(x(1))^2);

options = spot_sdp_default_options();
options.verbose = 1;
result = prog.minimize(0, solver,options);

Q=double(result.eval(Q))

V = result.eval(V);
Vdot = result.eval(Vdot)


[Theta,ThetaDot] = meshgrid(-pi:0.1:pi, -8:0.25:8);  % for surfs

figure(1);
subplot(1,4,1);
Vmesh = dmsubs(V,X,[sin(Theta(:)');cos(Theta(:)');ThetaDot(:)';sin(Theta(:)');cos(Theta(:)');ThetaDot(:)']);
surf(Theta,ThetaDot,reshape(Vmesh,size(Theta)));
title('$$ V $$','interpreter','latex','fontsize',20) 
xlabel('$$ \theta $$','interpreter','latex','fontsize',15)
ylabel('$$ \dot{\theta} $$','interpreter','latex','fontsize',15)

subplot(1,4,2);
Vmesh = dmsubs(V,X,[sin(Theta(:)');cos(Theta(:)');ThetaDot(:)';sin(Theta(:)');cos(Theta(:)');ThetaDot(:)']);
surf(sin(Theta),cos(Theta),reshape(Vmesh,size(Theta)));
title('$$ V $$','interpreter','latex','fontsize',20) 
xlabel('$$ \sin(\theta) $$','interpreter','latex','fontsize',15)
ylabel('$$ \cos(\theta) $$','interpreter','latex','fontsize',15)

subplot(1,4,3);
Vdotmesh = dmsubs(Vdot,X,[sin(Theta(:)');cos(Theta(:)');ThetaDot(:)';sin(Theta(:)');cos(Theta(:)');ThetaDot(:)']);
surf(Theta,ThetaDot,reshape(Vdotmesh,size(Theta)));
title('$$ Vdot $$','interpreter','latex','fontsize',20) 
xlabel('$$ \theta $$','interpreter','latex','fontsize',15)
ylabel('$$ \dot{\theta} $$','interpreter','latex','fontsize',15)

subplot(1,4,4);
Vdotmesh = dmsubs(Vdot,X,[sin(Theta(:)');cos(Theta(:)');ThetaDot(:)';sin(Theta(:)');cos(Theta(:)');ThetaDot(:)']);
surf(sin(Theta),cos(Theta),reshape(Vdotmesh,size(Theta)));
title('$$ Vdot $$','interpreter','latex','fontsize',20) 
xlabel('$$ \sin(\theta) $$','interpreter','latex','fontsize',15)
ylabel('$$ \cos(\theta) $$','interpreter','latex','fontsize',15)

% Vdotmesh = dmsubs(V,x,[sin(Theta(:)');cos(Theta(:)');ThetaDot(:)']);
% surf(Theta,ThetaDot,reshape(Vmesh,size(Theta)));
% title('$$ V $$','interpreter','latex','fontsize',20) 
% xlabel('$$ \theta $$','interpreter','latex','fontsize',15)
% ylabel('$$ \dot{\theta} $$','interpreter','latex','fontsize',15)

% figure(2);
% subplot(1,2,1);
% ezcontour(@(theta,thetadot)dmsubs(Vdot,x,[sin(theta');cos(theta');thetadot']),[-2*pi,2*pi,-8,8]);
% title('$$ \dot{V} $$','interpreter','latex','fontsize',20) 
% xlabel('$$ \theta $$','interpreter','latex','fontsize',15)
% ylabel('$$ \dot{\theta} $$','interpreter','latex','fontsize',15)
% subplot(1,2,2);
% Vdotmesh = dmsubs(Vdot,x,[sin(Theta(:)');cos(Theta(:)');ThetaDot(:)']);
% surf(Theta,ThetaDot,reshape(Vdotmesh,size(Theta)));
% title('$$ \dot{V} $$','interpreter','latex','fontsize',20) 
% xlabel('$$ \theta $$','interpreter','latex','fontsize',15)
% ylabel('$$ \dot{\theta} $$','interpreter','latex','fontsize',15)

end
