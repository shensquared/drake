function [V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,level,last_rho,delta_rho,last_V,all_V,do_plots)
	sum_rho=last_rho+delta_rho;
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);

	Vmonom = monomials(x,0:1);
	[prog,V] = prog.newFreePoly(Vmonom,4);
	
	Lmonom = monomials(x,0:2);
	[prog,L] = prog.newSOSPoly(Lmonom,16);

	vert0=[0;0];
	vert1=[0;last_rho];
	vert2=[-last_rho;0];
	vert3=[0;-last_rho];
	vert4=[last_rho;0];
	vert5=[0;last_rho+delta_rho];
	vert6=[-last_rho-delta_rho;0];
	vert7=[0;-last_rho-delta_rho];
	vert8=[last_rho+delta_rho;0];
	inner_verts=[vert1,vert2,vert3,vert4];
	outter_verts=[vert5,vert6,vert7,vert8];

	if level==0
		% all V zero at zero
		prog=prog.withEqs(subs(V,x,vert0)-zeros(4,1));
		% all V strictly positive at verts
		prog=prog.withPos((subs(V(1),x,vert5))-1e-6);
		prog=prog.withPos((subs(V(4),x,vert5))-1e-6);
		prog=prog.withPos((subs(V(1),x,vert6))-1e-6);
		prog=prog.withPos((subs(V(2),x,vert6))-1e-6);
		prog=prog.withPos((subs(V(2),x,vert7))-1e-6);
		prog=prog.withPos((subs(V(3),x,vert7))-1e-6);
		prog=prog.withPos((subs(V(3),x,vert8))-1e-6);
		prog=prog.withPos((subs(V(4),x,vert8))-1e-6);
	else
		last_max=max((dmsubs(last_V,x,inner_verts)),[],2);
		% evaluated at inner verts, last_V>=V
		prog=prog.withPos(last_max(1)-(subs(V(1),x,vert1)));
		prog=prog.withPos(last_max(1)-(subs(V(4),x,vert1)));
		prog=prog.withPos(last_max(2)-(subs(V(1),x,vert2)));
		prog=prog.withPos(last_max(2)-(subs(V(2),x,vert2)));
		prog=prog.withPos(last_max(3)-(subs(V(2),x,vert3)));
		prog=prog.withPos(last_max(3)-(subs(V(3),x,vert3)));
		prog=prog.withPos(last_max(4)-(subs(V(3),x,vert4)));
		prog=prog.withPos(last_max(4)-(subs(V(4),x,vert4)));
		% V at outter verts are bigger than inner verts
		prog=prog.withPos((subs(V(1),x,vert5-vert1)));
		prog=prog.withPos((subs(V(4),x,vert5-vert1)));
		prog=prog.withPos((subs(V(1),x,vert6-vert2)));
		prog=prog.withPos((subs(V(2),x,vert6-vert2)));
		prog=prog.withPos((subs(V(2),x,vert7-vert3)));
		prog=prog.withPos((subs(V(3),x,vert7-vert3)));
		prog=prog.withPos((subs(V(3),x,vert8-vert4)));
		prog=prog.withPos((subs(V(4),x,vert8-vert4)));
	end
	w=diff(V,x);
	% slack variables, pushing the solution into the interior of the feasible set
	[prog,slack]=prog.newPos(4);

	constraint1=[x(1);-x(2);-x(1)+x(2)-sum_rho;-(-x(1)+x(2)-last_rho)];
	constraint2=[x(1);x(2);-x(1)-x(2)-sum_rho;-(-x(1)-x(2)-last_rho)];
	constraint3=[x(2);-x(1);x(1)-x(2)-sum_rho;-(x(1)-x(2)-last_rho)];
	constraint4=[-x(1);-x(2);x(1)+x(2)-sum_rho;-(x(1)+x(2)-last_rho)];

	sum_rho_order=10^(floor(log10(sum_rho)));
	V1dot=w(1,:)*xdot*sum_rho_order;
	V2dot=w(2,:)*xdot*sum_rho_order;
	V3dot=w(3,:)*xdot*sum_rho_order;
	V4dot=w(4,:)*xdot*sum_rho_order;
	
	prog=prog.withSOS(-slack(1)-V1dot+L(1:4)'*constraint1);
	prog=prog.withSOS(-slack(2)-V2dot+L(5:8)'*constraint2);
	prog=prog.withSOS(-slack(3)-V3dot+L(9:12)'*constraint3);
	prog=prog.withSOS(-slack(4)-V4dot+L(13:16)'*constraint4);
	
	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(slack),@spot_mosek,options);

	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
        V=sol.eval(V);
        if level==0
        	all_V=V;
        else
        	all_V=[all_V;V];
        end
        plots_stuff(x,xdot,V,all_V,last_rho,delta_rho);        
        rho=sum_rho;
	else 
		sol_OK=false;
		% falls back to the last valid rho
		rho=last_rho;
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

