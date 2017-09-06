function VanDerPol_PWA()
checkDependency('spotless');
checkDependency('mosek');

x=msspoly('x',2);

%%% double cubic
% xdot = [-x(1)+x(1)^3; -x(2)+x(2)^3];

%%% VanDerPol
xdot = [x(2); -x(1)-x(2).*(x(1).^2-1)];
% df = [0 1; -1-2*x(1)*x(2), -(x(1)^2-1)];
A_zero=[0 1;-1 1];

[Vertices_values,w,old_rho]=diamond_iteration(x,xdot,A_zero);
plots(w);
rho=optimizeRho(x,xdot,Vertices_values,old_rho);
% w=level_set_version(x,xdot,A_zero);
% w=lp_version(x,xdot,A_zero);
end


function [Vertices_values,w,old_rho]=diamond_iteration(x,xdot,A_zero)
prog = spotsosprog;
prog = prog.withIndeterminate(x);
% Vertices conditions
[prog,Vertices_values] = prog.newPos(4);
% epsi=1e-4*ones(4,1);
% prog=prog.withEqs(Vertices_values(1)-1);

% prog=prog.withPos(Vertices_values(1)-1e-5);
rho=1e-4;

% the normals
w1=[-Vertices_values(2);Vertices_values(1)];
w2=[-Vertices_values(2);-Vertices_values(3)];
w3=[Vertices_values(4);-Vertices_values(3)];
w4=[Vertices_values(4);Vertices_values(1)];
% Lagrange multipliers
Lmonom = monomials(x,0:2);
[prog,L1] = prog.newFreePoly(Lmonom);
[prog,L2] = prog.newFreePoly(Lmonom);
[prog,L3] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L1);
prog = prog.withSOS(L2);
prog = prog.withSOS(L3);
[prog,L4] = prog.newFreePoly(Lmonom);
[prog,L5] = prog.newFreePoly(Lmonom);
[prog,L6] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L4);
prog = prog.withSOS(L5);
prog = prog.withSOS(L6);
[prog,L7] = prog.newFreePoly(Lmonom);
[prog,L8] = prog.newFreePoly(Lmonom);
[prog,L9] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L7);
prog = prog.withSOS(L8);
prog = prog.withSOS(L9);
[prog,L10] = prog.newFreePoly(Lmonom);
[prog,L11] = prog.newFreePoly(Lmonom);
[prog,L12] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L10);
prog = prog.withSOS(L11);
prog = prog.withSOS(L12);

% polytope 1
V1dot=w1'*xdot;
constraint1=x(1);
constraint2=-x(2);
constraint3=-x(1)+x(2)-rho;
prog=prog.withSOS(-V1dot+L1*constraint1+L2*constraint2+L3*constraint3);

% polytope 2
V2dot=w2'*xdot;
constraint4=x(1);
constraint5=x(2);
constraint6=-x(1)-x(2)-rho;
prog=prog.withSOS(-V2dot+L4*constraint4+L5*constraint5+L6*constraint6);

% polytope 3
V3dot=w3'*xdot;
constraint7=x(2);
constraint8=-x(1);
constraint9=x(1)-x(2)-rho;
prog=prog.withSOS(-V3dot+L7*constraint7+L8*constraint8+L9*constraint9);

% polytope 4
V4dot=w4'*xdot;
constraint10=-x(1);
constraint11=-x(2);
constraint12=x(1)+x(2)-rho;
prog=prog.withSOS(-V4dot+L10*constraint10+L11*constraint11+L12*constraint12);

options = spot_sdp_default_options();
options.verbose=0;
sol=prog.minimize(-sum(Vertices_values),@spot_mosek,options);
Vertices_values=double(sol.eval(Vertices_values))
w1=double((sol.eval(w1)));
w2=double(sol.eval(w2));
w3=double(sol.eval(w3));
w4=double(sol.eval(w4));
w=[w1,w2,w3,w4]
old_rho=min(Vertices_values);
end


function rho=optimizeRho(x,xdot,Vertices_values,old_rho)
prog = spotsosprog;
prog = prog.withIndeterminate(x);
% Vertices conditions
% [prog,rho] = prog.newPos(4);
% prog=prog.withPos(rho-old_rho*ones(4,1));
% [prog,infnorm] = prog.newPos(1);
% prog=prog.withPos(infnorm*ones(4,1)-rho);
% the normals
rho=old_rho*ones(4,1);
w1=[-Vertices_values(2);Vertices_values(1)];
w2=[-Vertices_values(2);-Vertices_values(3)];
w3=[Vertices_values(4);-Vertices_values(3)];
w4=[Vertices_values(4);Vertices_values(1)];

% Lagrange multipliers
Lmonom = monomials(x,0:2);
[prog,L1] = prog.newFreePoly(Lmonom);
[prog,L2] = prog.newFreePoly(Lmonom);
[prog,L3] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L1);
prog = prog.withSOS(L2);
prog = prog.withSOS(L3);
[prog,L4] = prog.newFreePoly(Lmonom);
[prog,L5] = prog.newFreePoly(Lmonom);
[prog,L6] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L4);
prog = prog.withSOS(L5);
prog = prog.withSOS(L6);
[prog,L7] = prog.newFreePoly(Lmonom);
[prog,L8] = prog.newFreePoly(Lmonom);
[prog,L9] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L7);
prog = prog.withSOS(L8);
prog = prog.withSOS(L9);
[prog,L10] = prog.newFreePoly(Lmonom);
[prog,L11] = prog.newFreePoly(Lmonom);
[prog,L12] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L10);
prog = prog.withSOS(L11);
prog = prog.withSOS(L12);

% polytope 1
V1dot=w1'*xdot;
constraint1=x(1);
constraint2=-x(2);
constraint3=w1'*x-rho(1);
prog=prog.withSOS(-L3*V1dot+L1*constraint1+L2*constraint2+constraint3);

% polytope 2
V2dot=w2'*xdot;
constraint4=x(1);
constraint5=x(2);
constraint6=w2'*x-rho(2);
prog=prog.withSOS(-L6*V2dot+L4*constraint4+L5*constraint5+constraint6);

% polytope 3
V3dot=w3'*xdot;
constraint7=x(2);
constraint8=-x(1);
constraint9=w3'*x-rho(3);
prog=prog.withSOS(-L9*V3dot+L7*constraint7+L8*constraint8+constraint9);

% polytope 4
V4dot=w4'*xdot;
constraint10=-x(1);
constraint11=-x(2);
constraint12=w4'*x-rho(4);
prog=prog.withSOS(-L12*V4dot+L10*constraint10+L11*constraint11+constraint12);

options = spot_sdp_default_options();
options.verbose=0;
sol=prog.minimize(-sum(rho),@spot_mosek,options);
rho=double(sol.eval(rho))
% w1=double((sol.eval(w1)));
% w2=double(sol.eval(w2));
% w3=double(sol.eval(w3));
% w4=double(sol.eval(w4));
% w=[w1,w2,w3,w4]

end


function square_iteration(x,xdot,delta,diamond_vert_values)
prog = spotsosprog;
prog = prog.withIndeterminate(x);
% Vertices conditions
[prog,Vertices_values] = prog.newPos(4);
% epsi=1e-4*ones(4,1);
% prog=prog.withPos(Vertices_values-epsi);
rho=1;

% the normals
w1=[-Vertices_values(2);Vertices_values(1)]./delta;
w2=[-Vertices_values(2);-Vertices_values(3)]./delta;
w3=[Vertices_values(4);-Vertices_values(3)]./delta;
w4=[Vertices_values(4);Vertices_values(1)]./delta;
% Lagrange multipliers
Lmonom = monomials(x,0:2);
[prog,L1] = prog.newFreePoly(Lmonom);
[prog,L2] = prog.newFreePoly(Lmonom);
[prog,L3] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L1);
prog = prog.withSOS(L2);
prog = prog.withSOS(L3);
[prog,L4] = prog.newFreePoly(Lmonom);
[prog,L5] = prog.newFreePoly(Lmonom);
[prog,L6] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L4);
prog = prog.withSOS(L5);
prog = prog.withSOS(L6);
[prog,L7] = prog.newFreePoly(Lmonom);
[prog,L8] = prog.newFreePoly(Lmonom);
[prog,L9] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L7);
prog = prog.withSOS(L8);
prog = prog.withSOS(L9);
[prog,L10] = prog.newFreePoly(Lmonom);
[prog,L11] = prog.newFreePoly(Lmonom);
[prog,L12] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L10);
prog = prog.withSOS(L11);
prog = prog.withSOS(L12);

% polytope 1
V1dot=w1'*xdot;
constraint1=x(1);
constraint2=-x(2);
constraint3=-x(1)+x(2)-rho;
prog=prog.withSOS(-V1dot+L1*constraint1+L2*constraint2+L3*constraint3);

% polytope 2
V2dot=w2'*xdot;
constraint4=x(1);
constraint5=x(2);
constraint6=-x(1)-x(2)-rho;
prog=prog.withSOS(-V2dot+L4*constraint4+L5*constraint5+L6*constraint6);

% polytope 3
V3dot=w3'*xdot;
constraint7=x(2);
constraint8=-x(1);
constraint9=x(1)-x(2)-rho;
prog=prog.withSOS(-V3dot+L7*constraint7+L8*constraint8+L9*constraint9);

% polytope 4
V4dot=w4'*xdot;
constraint10=-x(1);
constraint11=-x(2);
constraint12=x(1)+x(2)-rho;
prog=prog.withSOS(-V4dot+L10*constraint10+L11*constraint11+L12*constraint12);

options = spot_sdp_default_options();
options.verbose=0;
sol=prog.minimize(-sum(Vertices_values),@spot_mosek,options);
Vertices_values=double(sol.eval(Vertices_values));
w1=double((sol.eval(w1)));
w2=double(sol.eval(w2));
w3=double(sol.eval(w3));
w4=double(sol.eval(w4));
w=[w1,w2,w3,w4]
end

function plots(w,xdot)
x1=-1:.01:0;
y1=0:.01:1;
[X1,Y1]=meshgrid(x1,y1);
z1=w(1,1)*X1+w(2,1)*Y1;



x2=-1:.01:0;
y2=-1:.01:0;
[X2,Y2]=meshgrid(x2,y2);
z2=w(1,2)*X2+w(2,2)*Y2;

x3=0:.01:1;
y3=-1:.01:0;
[X3,Y3]=meshgrid(x3,y3);
z3=w(1,3)*X3+w(2,3)*Y3;

x4=0:.01:1;
y4=0:.01:1;
[X4,Y4]=meshgrid(x4,y4);
z4=w(1,4)*X4+w(2,4)*Y4;


[fullX,fullY]=meshgrid(-2:.05:2,-2:.05:2);
subplot(1,2,1);
surf(X1,Y1,z1);hold on
surf(X2,Y2,z2); hold on
surf(X3,Y3,z3);hold on
surf(X4,Y4,z4);

subplot(1,2,2);
contour(X1,Y1,z1);hold on
contour(X2,Y2,z2);hold on
contour(X3,Y3,z3);hold on
contour(X4,Y4,z4);hold on
quiver(fullX,fullY,fullY,-fullX-fullY.*(fullX.^2-1))

% subplot(1,3,3);
% surf(X1,Y1,z1);hold on
% surf(X2,Y2,z2);hold on
% surf(X3,Y3,z3);hold on
% surf(X4,Y4,z4);hold on


% figure(2);
% tri = delaunay(X1,Y1);
% trisurf(tri,X1,Y1,z1);
% X=[X1,X2,X3,X4];
% Y=[Y1,Y2,Y3,Y4];
% z=[z1,z2,z3,z4];
% surf(X,Y,z);
% contour(X,Y,z);
end

function dot=get_dot(x)
	dot=[x(2); -x(1)-x(2).*(x(1).^2-1)];
end
