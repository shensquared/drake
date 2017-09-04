function Vertices_values=VanDerPol_PWA()
checkDependency('spotless');
checkDependency('mosek');
x=msspoly('x',2);

%%% double cubic
% xdot = [-x(1)+x(1)^3; -x(2)+x(2)^3];

%%% VanDerPol
xdot = [x(2); -x(1)-x(2).*(x(1).^2-1)];
% df = [0 1; -1-2*x(1)*x(2), -(x(1)^2-1)];
A_zero=[0 1;-1 1];

% epsi=1e-4*ones(4,1);
w=sos_version(x,xdot,A_zero);
% w=lp_version(x,xdot,A_zero);
plots(w);
end


function w=sos_version(x,xdot,A_zero)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	% Vertices conditions
	[prog,Vertices_values] = prog.newPos(4);

	% prog=prog.withPos(Vertices_values-epsi);
	% the normals
	w1=[-Vertices_values(2);Vertices_values(1)];
	w2=[-Vertices_values(2);-Vertices_values(3)];
	w3=[Vertices_values(4);-Vertices_values(3)];
	w4=[Vertices_values(4);Vertices_values(1)];

	Lmonom = monomials(x,0:2);

	[prog,L1] = prog.newFreePoly(Lmonom);
	[prog,L2] = prog.newFreePoly(Lmonom);
	[prog,L3] = prog.newFreePoly(Lmonom);

	prog = prog.withSOS(L1);
	prog = prog.withSOS(L2);
	prog = prog.withSOS(L3);

	% polytope 1
	V1dot=w1'*xdot;
	constraint1=x(1);
	constraint2=-x(2);
	constraint3=-x(1)+x(2)-1e-0;
	prog=prog.withSOS(-V1dot+L1*constraint1+L2*constraint2+L3*constraint3);

	[prog,L1] = prog.newFreePoly(Lmonom);
	[prog,L2] = prog.newFreePoly(Lmonom);
	[prog,L3] = prog.newFreePoly(Lmonom);

	prog = prog.withSOS(L1);
	prog = prog.withSOS(L2);
	prog = prog.withSOS(L3);

	V1dot=w1'*xdot;
	constraint1=x(1);
	constraint2=-x(2);
	constraint3=-x(1)+x(2);
	prog=prog.withSOS(-V1dot+L1*constraint1+L2*constraint2+L3*constraint3);


	% polytope 2
	[prog,L4] = prog.newFreePoly(Lmonom);
	[prog,L5] = prog.newFreePoly(Lmonom);
	[prog,L6] = prog.newFreePoly(Lmonom);

	prog = prog.withSOS(L4);
	prog = prog.withSOS(L5);
	prog = prog.withSOS(L6);

	V2dot=w2'*xdot;
	constraint4=x(1);
	constraint5=x(2);
	constraint6=-x(1)-x(2)-1e-0;
	prog=prog.withSOS(-V2dot+L4*constraint4+L5*constraint5+L6*constraint6);


	% polytope 3
	[prog,L7] = prog.newFreePoly(Lmonom);
	[prog,L8] = prog.newFreePoly(Lmonom);
	[prog,L9] = prog.newFreePoly(Lmonom);

	prog = prog.withSOS(L7);
	prog = prog.withSOS(L8);
	prog = prog.withSOS(L9);

	V3dot=w3'*xdot;
	constraint7=x(2);
	constraint8=-x(1);
	constraint9=x(1)-x(2)-1e-0;
	prog=prog.withSOS(-V3dot+L7*constraint7+L8*constraint8+L9*constraint9);

	% polytope 4
	[prog,L10] = prog.newFreePoly(Lmonom);
	[prog,L11] = prog.newFreePoly(Lmonom);
	[prog,L12] = prog.newFreePoly(Lmonom);
	prog = prog.withSOS(L10);
	prog = prog.withSOS(L11);
	prog = prog.withSOS(L12);

	V4dot=w4'*xdot;
	constraint10=-x(1);
	constraint11=-x(2);
	constraint12=x(1)+x(2)-1e-0;
	prog=prog.withSOS(-V4dot+L10*constraint10+L11*constraint11+L12*constraint12);


	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(Vertices_values(2)),@spot_mosek,options);
	Vertices_values=double(sol.eval(Vertices_values));
	w1=double((sol.eval(w1)));
	w2=double(sol.eval(w2));
	w3=double(sol.eval(w3));
	w4=double(sol.eval(w4));
	w=[w1,w2,w3,w4]
end

function plots(w)
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

	% X=[X1,X2,X3,X4];
	% Y=[Y1,Y2,Y3,Y4];
	% z=[z1,z2,z3,z4];
	% surf(X,Y,z);
	% contour(X,Y,z);
end


function w=lp_version(x,xdot,A_zero);
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	% Vertices conditions
	[prog,Vertices_values] = prog.newPos(4);
	% prog=prog.withPos(Vertices_values-epsi);
	
	% the normals
	w1=[-Vertices_values(2);Vertices_values(1)];
	w2=[-Vertices_values(2);-Vertices_values(3)];
	w3=[Vertices_values(4);-Vertices_values(3)];
	w4=[Vertices_values(4);Vertices_values(1)];

	prog = prog.withPos(-w1'*A_zero*([0;1]));
	prog = prog.withPos(-w1'*A_zero*([-1;0]));

	prog = prog.withPos(-w2'*A_zero*([-1;0]));
	prog = prog.withPos(-w2'*A_zero*([0;-1]));

	prog = prog.withPos(-w3'*A_zero*([0;-1]));
	prog = prog.withPos(-w3'*A_zero*([1;0]));

	prog = prog.withPos(-w4'*A_zero*([1;0]));
	prog = prog.withPos(-w4'*A_zero*([0;1]));

	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(Vertices_values),@spot_mosek,options);
	Vertices_values=double(sol.eval(Vertices_values))
	w1=double((sol.eval(w1)));
	w2=double(sol.eval(w2));
	w3=double(sol.eval(w3));
	w4=double(sol.eval(w4));
	w=[w1,w2,w3,w4]
end


function rho=find_boundary(V,x,options,old_rho)
prog = spotsosprog;
prog = prog.withIndeterminate(x);
Vmonom = monomials(x,0:1);
[prog,V] = prog.newFreePoly(Vmonom);

Lmonom = monomials(x,0:2);
[prog,L1] = prog.newFreePoly(Lmonom);
[prog,L2] = pog.newFreePoly(Lmonom);

[prog,rho] = prog.newPos(1);
% prog=prog.withPos(rho-old_rho);
Vdot=diff(V,x)*xdot;
constraint1=x(1)*(rho-x(1));
constraint2=x(2)*(rho-x(2));
prog=prog.withSOS(-Vdot+L1*constraint1+L2*constraint2);

options = spot_sdp_default_options();
options.verbose=0;
sol=prog.minimize(-rho,@spot_mosek,options);
V=sol.eval(V);
end