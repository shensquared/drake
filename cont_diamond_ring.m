function [rho,Vertices_values,w,sol_OK]=cont_diamond_ring(x,xdot,last_rho,last_vertice_values,delta_rho,do_plots)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	[prog,this_vert_values] = prog.newPos(4);
	prog=prog.withPos(this_vert_values-1e-6*last_rho*ones(4,1));
	Vmonom = monomials(x,0:1);
	[prog,V] = prog.newFreePoly(Vmonom,4);

	vert1=[0;last_rho];
	vert2=[-last_rho;0];
	vert3=[0;-last_rho];
	vert4=[last_rho;0];
	vert5=[0;last_rho+delta_rho];
	vert6=[-last_rho-delta_rho;0];
	vert7=[0;-last_rho-delta_rho];
	vert8=[last_rho+delta_rho;0];

	prog=prog.withPos(-(subs(V(1),x,vert1)-last_vertice_values(1)));
	prog=prog.withPos(-(subs(V(1),x,vert2)-last_vertice_values(2)));
	prog=prog.withPos(-(subs(V(2),x,vert2)-last_vertice_values(2)));
	prog=prog.withPos(-(subs(V(2),x,vert3)-last_vertice_values(3)));
	prog=prog.withPos(-(subs(V(3),x,vert3)-last_vertice_values(3)));
	prog=prog.withPos(-(subs(V(3),x,vert4)-last_vertice_values(4)));
	prog=prog.withPos(-(subs(V(4),x,vert4)-last_vertice_values(4)));
	prog=prog.withPos(-(subs(V(4),x,vert1)-last_vertice_values(1)));


	prog=prog.withPos((subs(V(1),x,vert5)-last_vertice_values(1)));
	prog=prog.withPos((subs(V(1),x,vert6)-last_vertice_values(2)));
	prog=prog.withPos((subs(V(2),x,vert6)-last_vertice_values(2)));
	prog=prog.withPos((subs(V(2),x,vert7)-last_vertice_values(3)));
	prog=prog.withPos((subs(V(3),x,vert7)-last_vertice_values(3)));
	prog=prog.withPos((subs(V(3),x,vert8)-last_vertice_values(4)));
	prog=prog.withPos((subs(V(4),x,vert8)-last_vertice_values(4)));
	prog=prog.withPos((subs(V(4),x,vert5)-last_vertice_values(1)));


	% prog=prog.withEqs(subs(V(1),x,vert5)-this_vert_values(1));
	% prog=prog.withEqs(subs(V(1),x,vert5)-this_vert_values(1));
	% prog=prog.withEqs(subs(V(1),x,vert6)-this_vert_values(2));
	% prog=prog.withEqs(subs(V(2),x,vert6)-this_vert_values(2));
	% prog=prog.withEqs(subs(V(2),x,vert7)-this_vert_values(3));
	% prog=prog.withEqs(subs(V(3),x,vert7)-this_vert_values(3));
	% prog=prog.withEqs(subs(V(3),x,vert8)-this_vert_values(4));
	% prog=prog.withEqs(subs(V(4),x,vert8)-this_vert_values(4));
	% prog=prog.withEqs(subs(V(4),x,vert5)-this_vert_values(1));
	
	Lmonom = monomials(x,0:2);
	[prog,L] = prog.newSOSPoly(Lmonom,16);

	sum_rho=last_rho+delta_rho;
	r_sum_rho=1/sum_rho;
	r_delta_rho=1/delta_rho;
	r_last_rho=1/last_rho;

	% % slack variables, pushing the solution into the interior of the feasible set
	[prog,slack]=prog.newPos(4);

	constraint1=[x(1);-x(2);-x(1)+x(2)-sum_rho;-(-x(1)+x(2)-last_rho)];
	constraint2=[x(1);x(2);-x(1)-x(2)-sum_rho;-(-x(1)-x(2)-last_rho)];
	constraint3=[x(2);-x(1);x(1)-x(2)-sum_rho;-(x(1)-x(2)-last_rho)];
	constraint4=[-x(1);-x(2);x(1)+x(2)-sum_rho;-(x(1)+x(2)-last_rho)];
	% manufactured scaling, somehow improves numerical condition...
	V1dot=diff(V(1),x)*xdot*sum_rho;
	V2dot=diff(V(2),x)*xdot*sum_rho;
	V3dot=diff(V(3),x)*xdot*sum_rho;
	V4dot=diff(V(4),x)*xdot*sum_rho;
	prog=prog.withSOS(-slack(1)-V1dot+L(1:4)'*constraint1);
	prog=prog.withSOS(-slack(2)-V2dot+L(5:8)'*constraint2);
	prog=prog.withSOS(-slack(3)-V3dot+L(9:12)'*constraint3);
	prog=prog.withSOS(-slack(4)-V4dot+L(13:16)'*constraint4);

	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(0),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		w1=double(sol.eval(diff(V(1),x)));
		w2=double(sol.eval(diff(V(2),x)));
		w3=double(sol.eval(diff(V(3),x)));
		w4=double(sol.eval(diff(V(4),x)));
		w=[w1,w2,w3,w4];
		delta_rho=double(sol.eval(delta_rho));
		b1=double(sol.eval(subs(V(1),x,zeros(2,1))));
		b2=double(sol.eval(subs(V(2),x,zeros(2,1))));
		b3=double(sol.eval(subs(V(3),x,zeros(2,1))));
		b4=double(sol.eval(subs(V(4),x,zeros(2,1))));
		% b1=0;
		% b2=0;
		% b3=0;
		% b4=0;

		V1=(sol.eval(V(1)));
		V2=(sol.eval(V(2)));
		V3=(sol.eval(V(3)));
		V4=(sol.eval(V(4)));

		Vertices_values=double(sol.eval(this_vert_values))
        sol.eval(V1dot)
        rho=delta_rho+last_rho;
		% rho=double(sol.eval(rho));
		% w=[w1,w2,w3,w4]./rho;
	else
		sol_OK=false;
	end

	if do_plots
		% figure(2)
		syms x1 x2;
		V= piecewise(x1<=0&-x2<=0&-x1+x2-sum_rho<=0&(-(-x1+x2-last_rho))<=0,w1*[x1;x2]+b1,x1<=0&x2<=0&-x1-x2-sum_rho<=0&(-(-x1-x2-last_rho))<=0,w2*[x1;x2]+b2,x2<=0&-x1<=0&x1-x2-sum_rho<=0&(-(x1-x2-last_rho))<=0,w3*[x1;x2]+b3,-x1<=0&-x2<=0&x1+x2-sum_rho<=0&(-(x1+x2-last_rho))<=0,w4*[x1;x2]+b4,NaN);
		% V= piecewise(x1<=0&-x2<=0&-x1+x2-rho<=0&(-(-x1+x2-last_rho))<=0,10,x1<=0&x2<=0&-x1-x2-rho<=0&(-(-x1-x2-last_rho))<=0,10,x2<=0&-x1<=0&x1-x2-rho<=0&(-(x1-x2-last_rho))<=0,10,-x1<=0&-x2<=0&x1+x2-rho<=0&(-(x1+x2-last_rho))<=0,10,NaN);
		% subplot(1,2,1)
		fsurf(V); hold on
		disp('finished plotting');
		% subplot(1,2,2)
		% b1=0;
		% b2=0;
		% b3=0;
		% b4=0;
		% V= piecewise(x1<=0&-x2<=0&-x1+x2-sum_rho<=0&(-(-x1+x2-last_rho))<=0,w1*[x1;x2]+b1,x1<=0&x2<=0&-x1-x2-sum_rho<=0&(-(-x1-x2-last_rho))<=0,w2*[x1;x2]+b2,x2<=0&-x1<=0&x1-x2-sum_rho<=0&(-(x1-x2-last_rho))<=0,w3*[x1;x2]+b3,-x1<=0&-x2<=0&x1+x2-sum_rho<=0&(-(x1+x2-last_rho))<=0,w4*[x1;x2]+b4,NaN);
		% fsurf(V);
		% fcontour(V);
		% figure(3)
		% % sym_xdot=[-2*x1+x1^3; -2*x2+x2^3];
		% sym_xdot =-[x2; -x1-x2.*(x1.^2-1)];
		% Vdot=piecewise(x1<=0&-x2<=0&-x1+x2-rho<=0,w(:,1)'*sym_xdot,x1<=0&x2<=0&-x1-x2-rho<=0,w(:,2)'*sym_xdot,x2<=0&-x1<=0&x1-x2-rho<=0,w(:,3)'*sym_xdot,-x1<=0&-x2<=0&x1+x2-rho<=0,w(:,4)'*sym_xdot,NaN);
		% fsurf(Vdot); hold on
	end

end
