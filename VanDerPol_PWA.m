function VanDerPol_PWA()
	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);
	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	xdot = [x(2); -x(1)-x(2).*(x(1).^2-1)];
	df = [0 1; -1-2*x(1)*x(2), -(x(1)^2-1)];

	rhomin=0; rhomax=100;
	rho=2e-3;
	iter_count=0;
	while(iter_count<100)
	% 	 disp(iter_count);
		[Vertices_values,w,rho,L]=diamond(x,xdot,df,'fix_rho',rho);
		rho=1.2*rho;
		 % [Vertices_values,w,rho,L,OK_flag]=diamond(x,xdot,df,'fix_V_L',Vertices_values,L,w,rho);
		iter_count=iter_count+1;
	end
	% plots(w);
end

function [Vertices_values,w,rho,L]=diamond(x,xdot,df,method,varargin)
	% disp(method);
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	switch method
	case 'fix_rho'
		rho=varargin{1};
		% xdot=rho*subs(xdot,x,1/rho*x);
		[prog,Vertices_values] = prog.newPos(4);
		prog=prog.withPos(Vertices_values()-1e-6*ones(4,1));
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
		% the normals, without scaling (so that w is at the same scale as vertice values)
		w1=[-Vertices_values(2);Vertices_values(1)];
		w2=[-Vertices_values(2);-Vertices_values(3)];
		w3=[Vertices_values(4);-Vertices_values(3)];
		w4=[Vertices_values(4);Vertices_values(1)];
		% slack variables, pushing the solution into the interior of the feasible set
		[prog,slack]=prog.newPos(4);

	case 'fix_V_L'
		Vertices_values=varargin{1};
		L=varargin{2};
		w=varargin{3};
		old_rho=varargin{4};
		L1=L(1);
		L2=L(2);
		L3=L(3);
		L4=L(4);
		L5=L(5);
		L6=L(6);
		L7=L(7);
		L8=L(8);
		L9=L(9);
		L10=L(10);
		L11=L(11);
		L12=L(12);
		% the decision variable is rho
		[prog,rho] = prog.newPos(1);
		prog=prog.withPos(rho-old_rho);
		w1=w(:,1);
		w2=w(:,2);
		w3=w(:,3);
		w4=w(:,4);
		[prog,slack]=prog.newPos(4);
	otherwise
		error('unknown method');
	end

	% polytope 1
	V1dot=w1'*xdot;
	constraint1=x(1);
	constraint2=-x(2);
	constraint3=-x(1)+x(2)-rho;
	prog=prog.withSOS(-slack(1)-V1dot+L1*constraint1+L2*constraint2+L3*constraint3);

	% polytope 2
	V2dot=w2'*xdot;
	constraint4=x(1);
	constraint5=x(2);
	constraint6=-x(1)-x(2)-rho;
	prog=prog.withSOS(-slack(2)-V2dot+L4*constraint4+L5*constraint5+L6*constraint6);

	% polytope 3
	V3dot=w3'*xdot;
	constraint7=x(2);
	constraint8=-x(1);
	constraint9=x(1)-x(2)-rho;
	prog=prog.withSOS(-slack(3)-V3dot+L7*constraint7+L8*constraint8+L9*constraint9);

	% polytope 4
	V4dot=w4'*xdot;
	constraint10=-x(1);
	constraint11=-x(2);
	constraint12=x(1)+x(2)-rho;
	prog=prog.withSOS(-slack(4)-V4dot+L10*constraint10+L11*constraint11+L12*constraint12);
	options = spot_sdp_default_options();
	options.verbose=0;

	% vert0=[0;0];
	% vert1=[0;rho];
	% vert2=[-rho;0];
	% vert3=[0;-rho];
	% vert4=[rho;0];
	% v1dot_at_verticies=[subs(V1dot,x,vert1);subs(V1dot,x,vert2);subs(V2dot,x,vert2);subs(V2dot,x,vert3);subs(V3dot,x,vert3);subs(V3dot,x,vert4);subs(V4dot,x,vert4);subs(V4dot,x,vert1)];

	sol=prog.minimize(-sum(slack),@spot_mosek,options);
	if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
	  error('The solver threw an internal error.');
	end
	if ~sol.isPrimalFeasible
		error('Problem looks primal infeasible.');
	end

	if ~sol.isDualFeasible
		error('Problem looks dual infeasible. It is probably unbounded.');
	end

	if sol.status~=spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		error('not primal and dual feasible');
	end

	Vertices_values=double(sol.eval(Vertices_values))
	w1=double((sol.eval(w1)));
	w2=double(sol.eval(w2));
	w3=double(sol.eval(w3));
	w4=double(sol.eval(w4));
	w=[w1,w2,w3,w4]
	rho=double(sol.eval(rho))
	% L=[(L1),(L2),(L3),(L4),(L5),(L6),(L7),(L8),(L9),(L10),(L11),(L12)];
	L=[sol.eval(L1),sol.eval(L2),sol.eval(L3),sol.eval(L4),sol.eval(L5),sol.eval(L6),sol.eval(L7),sol.eval(L8),sol.eval(L9),sol.eval(L10),sol.eval(L11),sol.eval(L12)];
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
