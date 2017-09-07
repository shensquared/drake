function square_iteration(x,xdot,diamond_vert_values)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	% Vertices conditions
	[prog,Vertices_values] = prog.newPos(4);
	% epsi=1e-4*ones(4,1);
	prog=prog.withPos(Vertices_values(1)-diamond_vert_values(1));
	prog=prog.withPos(Vertices_values(1)-diamond_vert_values(2));
	prog=prog.withPos(Vertices_values(2)-diamond_vert_values(2));
	prog=prog.withPos(Vertices_values(2)-diamond_vert_values(3));
	prog=prog.withPos(Vertices_values(3)-diamond_vert_values(3));
	prog=prog.withPos(Vertices_values(3)-diamond_vert_values(4));
	prog=prog.withPos(Vertices_values(4)-diamond_vert_values(4));
	prog=prog.withPos(Vertices_values(4)-diamond_vert_values(1));

	rho=1e-4;

	% the normals
	w1=[diamond_vert_values(1)-Vertices_values(1);Vertices_values(1)-diamond_vert_values(2)];
	w2=[diamond_vert_values(3)-Vertices_values(2);diamond_vert_values(2)-Vertices_values(2)];
	w3=[Vertices_values(3)-diamond_vert_values(3);diamond_vert_values(4)-Vertices_values(3)];
	w4=[Vertices_values(4)-diamond_vert_values(1);Vertices_values(4)-diamond_vert_values(4)];

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
	constraint1=-rho-x(1);
	constraint2=x(2)-rho;
	constraint3=-(-x(1)+x(2)-rho);
	prog=prog.withSOS(-V1dot+L1*constraint1+L2*constraint2+L3*constraint3);

	% polytope 2
	V2dot=w2'*xdot;
	constraint4=-rho-x(1);
	constraint5=-rho-x(2);
	constraint6=-(-x(1)-x(2)-rho);
	prog=prog.withSOS(-V2dot+L4*constraint4+L5*constraint5+L6*constraint6);

	% polytope 3
	V3dot=w3'*xdot;
	constraint7=-rho-x(2);
	constraint8=x(1)-rho;
	constraint9=-(x(1)-x(2)-rho);
	prog=prog.withSOS(-V3dot+L7*constraint7+L8*constraint8+L9*constraint9);

	% polytope 4
	V4dot=w4'*xdot;
	constraint10=-rho+x(1);
	constraint11=x(2)-rho;
	constraint12=-(x(1)+x(2)-rho);
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
function w=level_set_version(x,xdot,A_zero)
prog = spotsosprog;
prog = prog.withIndeterminate(x);
% Vertices conditions
[prog,Vertices_values] = prog.newPos(4);
% epsi=1e-4*ones(4,1);
% prog=prog.withPos(Vertices_values(1)-1e-5);
rho=1;

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
[prog,L4] = prog.newFreePoly(Lmonom);
[prog,L5] = prog.newFreePoly(Lmonom);
[prog,L6] = prog.newFreePoly(Lmonom);
[prog,L7] = prog.newFreePoly(Lmonom);
[prog,L8] = prog.newFreePoly(Lmonom);
[prog,L9] = prog.newFreePoly(Lmonom);
[prog,L10] = prog.newFreePoly(Lmonom);
[prog,L11] = prog.newFreePoly(Lmonom);
[prog,L12] = prog.newFreePoly(Lmonom);

prog = prog.withSOS(L1);
prog = prog.withSOS(L2);
prog = prog.withSOS(L3);
prog = prog.withSOS(L4);
prog = prog.withSOS(L5);
prog = prog.withSOS(L6);
prog = prog.withSOS(L7);
prog = prog.withSOS(L8);
prog = prog.withSOS(L9);
prog = prog.withSOS(L10);
prog = prog.withSOS(L11);
prog = prog.withSOS(L12);

% polytope 1
V1dot=w1'*xdot;
constraint1=w1'*x-1;
% constraint2=-x(2);
% constraint3=-x(1)+x(2)-rho;
prog=prog.withSOS(-V1dot+L1*constraint1);

% polytope 2
V2dot=w2'*xdot;
constraint4=w2'*x-1;
% constraint5=x(2);
% constraint6=-x(1)-x(2)-rho;
prog=prog.withSOS(-V2dot+L4*constraint4);

% polytope 3
V3dot=w3'*xdot;
constraint7=w3'*x-1;
% constraint8=-x(1);
% constraint9=x(1)-x(2)-rho;
prog=prog.withSOS(-V3dot+L7*constraint7);

% polytope 4
V4dot=w4'*xdot;
constraint10=w4'*x-1;
% constraint11=-x(2);
% constraint12=x(1)+x(2)-rho;
prog=prog.withSOS(-V4dot+L10*constraint10);

options = spot_sdp_default_options();
options.verbose=0;
sol=prog.minimize(-sum(Vertices_values),@spot_mosek,options);
Vertices_values=double(sol.eval(Vertices_values))
w1=double((sol.eval(w1)));
w2=double(sol.eval(w2));
w3=double(sol.eval(w3));
w4=double(sol.eval(w4));
w=[w1,w2,w3,w4]
end
