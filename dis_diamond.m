function [V,rho]=dis_diamond(x,xdot,do_plots,rho)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);

	Vmonom = monomials(x,0:1);
	[prog,V] = prog.newFreePoly(Vmonom,4);

	Lmonom = monomials(x,0:2);
	[prog,L] = prog.newSOSPoly(Lmonom,12);

	vert0=[0;0];
	vert1=[0;rho];
	vert2=[-rho;0];
	vert3=[-rho;0];
	vert4=[0;-rho];
	vert5=[0;-rho];
	vert6=[rho;0];
	vert7=[rho;0];
	vert8=[0;rho];

	prog=prog.withEqs(subs(V,x,vert0)-zeros(4,1));
% all V strictly positive at verts
	prog=prog.withPos((subs(V(1),x,vert1))-1e-6);
	prog=prog.withPos((subs(V(1),x,vert2))-1e-6);
	prog=prog.withPos((subs(V(2),x,vert3))-1e-6);
	prog=prog.withPos((subs(V(2),x,vert4))-1e-6);
	prog=prog.withPos((subs(V(3),x,vert5))-1e-6);
	prog=prog.withPos((subs(V(3),x,vert6))-1e-6);
	prog=prog.withPos((subs(V(4),x,vert7))-1e-6);
	prog=prog.withPos((subs(V(4),x,vert8))-1e-6);

% the normals, without scaling (so that w is at the same scale as vertice values)
	w=diff(V,x);
% slack variables, pushing the solution into the interior of the feasible set
	[prog,slack]=prog.newPos(4);

	constraint1=[x(1);-x(2);-x(1)+x(2)-rho;];
	constraint2=[x(1);x(2);-x(1)-x(2)-rho;];
	constraint3=[x(2);-x(1);x(1)-x(2)-rho;];
	constraint4=[-x(1);-x(2);x(1)+x(2)-rho;];

	V1dot=w(1,:)*xdot*rho;
	V2dot=w(2,:)*xdot*rho;
	V3dot=w(3,:)*xdot*rho;
	V4dot=w(4,:)*xdot*rho;

	prog=prog.withSOS(-slack(1)-V1dot+L(1:3)'*constraint1);
	prog=prog.withSOS(-slack(2)-V2dot+L(4:6)'*constraint2);
	prog=prog.withSOS(-slack(3)-V3dot+L(7:9)'*constraint3);
	prog=prog.withSOS(-slack(4)-V4dot+L(10:12)'*constraint4);

	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(slack),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
        V=sol.eval(V);
        w=sol.eval(w);
	else 
		sol_OK=false;
		% falls back to the last valid rho
		if do_plots
			figure(2)
			syms x1 x2;
			V= piecewise(x1<=0&-x2<=0&-x1+x2-rho<=0,w(:,1)'*[x1;x2],x1<=0&x2<=0&-x1-x2-rho<=0,w(:,2)'*[x1;x2],x2<=0&-x1<=0&x1-x2-rho<=0,w(:,3)'*[x1;x2],-x1<=0&-x2<=0&x1+x2-rho<=0,w(:,4)'*[x1;x2],NaN);
			fsurf(V); hold on
			% fcontour(V);
			figure(3)
			% sym_xdot=[-2*x1+x1^3; -2*x2+x2^3];
			sym_xdot =-[x2; -x1-x2.*(x1.^2-1)];
			Vdot=piecewise(x1<=0&-x2<=0&-x1+x2-rho<=0,w(:,1)'*sym_xdot,x1<=0&x2<=0&-x1-x2-rho<=0,w(:,2)'*sym_xdot,x2<=0&-x1<=0&x1-x2-rho<=0,w(:,3)'*sym_xdot,-x1<=0&-x2<=0&x1+x2-rho<=0,w(:,4)'*sym_xdot,NaN);
			fsurf(Vdot); hold on
		end
	end
end
