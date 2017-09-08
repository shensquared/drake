function VanDerPol_PWA()
	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);
	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	xdot = [x(2); -x(1)-x(2).*(x(1).^2-1)];
	df = [0 1; -1-2*x(1)*x(2), -(x(1)^2-1)];

	rho=5e-3;
	iter_count=0;
	Vertices_values=zeros(4,1);
	w=zeros(2,4);
	sol_OK=true;
	max_iter=100;
	epsi=1e-6;

	last_vertice_values =   1.0e-05 *[
    0.1232;
    0.1576;
    0.1227;
    0.1590];

	while(sol_OK)
		disp(iter_count);
		[rho,Vertices_values,w,sol_OK]=diamond(x,xdot,df,'fix_rho',rho,Vertices_values,w)
		rho=1.2*rho;
		iter_count=iter_count+1;
    end

    % beginning of the diamond ring iteration
	% while(sol_OK)
	% 	disp(iter_count);
% 	diamond_ring(x,xdot,df,0,zeros(4,1),5e-3);

		% rho=1.2*rho;
		% iter_count=iter_count+1;
    % end
	rho;
	% plots(w);
end

function diamond_ring(x,xdot,df,last_rho,last_vertice_values,delta_rho)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);


	[prog,t] = prog.newPos(1);
	prog=prog.withPos(t*ones(4,1)-last_vertice_values);
	[prog,this_vert_values] = prog.newPos(4);
	prog=prog.withPos(this_vert_values-last_vertice_values);

	Vmonom = monomials(x,0:1);
	[prog,V1] = prog.newFreePoly(Vmonom);
	[prog,V2] = prog.newFreePoly(Vmonom);
	[prog,V3] = prog.newFreePoly(Vmonom);
	[prog,V4] = prog.newFreePoly(Vmonom);

	vert1=[0;last_rho];
	vert2=[-last_rho;0];
	vert3=[0;-last_rho];
	vert4=[last_rho;0];
	vert5=[0;last_rho+delta_rho];
	vert6=[-last_rho-delta_rho;0];
	vert7=[0;-last_rho-delta_rho];
	vert8=[last_rho+delta_rho;0];

	prog=prog.withPos(subs(V1,x,vert1)-last_vertice_values(1));
	prog=prog.withPos(subs(V1,x,vert2)-last_vertice_values(2));
	prog=prog.withPos(subs(V2,x,vert2)-last_vertice_values(2));
	prog=prog.withPos(subs(V2,x,vert3)-last_vertice_values(3));
	prog=prog.withPos(subs(V3,x,vert3)-last_vertice_values(3));
	prog=prog.withPos(subs(V3,x,vert4)-last_vertice_values(4));
	prog=prog.withPos(subs(V4,x,vert4)-last_vertice_values(4));
	prog=prog.withPos(subs(V4,x,vert1)-last_vertice_values(1));

	% prog=prog.withEqs(subs(V1,x,vert5)-this_vert_values(1));

	prog=prog.withEqs(subs(V1,x,vert5)-this_vert_values(1));
	prog=prog.withEqs(subs(V1,x,vert6)-this_vert_values(2));
	prog=prog.withEqs(subs(V2,x,vert6)-this_vert_values(2));
	prog=prog.withEqs(subs(V2,x,vert7)-this_vert_values(3));
	prog=prog.withEqs(subs(V3,x,vert7)-this_vert_values(3));
	prog=prog.withEqs(subs(V3,x,vert8)-this_vert_values(4));
	prog=prog.withEqs(subs(V4,x,vert8)-this_vert_values(4));
	prog=prog.withEqs(subs(V4,x,vert5)-this_vert_values(1));
	
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
	% ring 1
	V1dot=diff(V1,x)*xdot*sum_rho;
	% constraint1=x(1);
	% constraint2=-x(2);
	% constraint3=-x(1)+x(2)-sum_rho;
	% constraint4=-(-x(1)+x(2)-last_rho);	
	prog=prog.withSOS(-slack(1)-V1dot+L(1:4)*constraint1);

	% ring 2
	V2dot=diff(V2,x)*xdot*sum_rho;
	% constraint4=x(1);
	% constraint5=x(2);
	% constraint6=-x(1)-x(2)-sum_rho;
	prog=prog.withSOS(-slack(2)-V2dot+L(5:8)*constraint2);

	% ring 3
	V3dot=diff(V3,x)*xdot*sum_rho;
	% constraint7=x(2);
	% constraint8=-x(1);
	% constraint9=x(1)-x(2)-sum_rho;
	prog=prog.withSOS(-slack(3)-V3dot+L(9:12)*constraint3);

	% ring 4
	V4dot=diff(V4,x)*xdot*sum_rho;
	% constraint10=-x(1);
	% constraint11=-x(2);
	% constraint12=x(1)+x(2)-sum_rho;
	prog=prog.withSOS(-slack(4)-V4dot+L(13:16)*constraint4);
	

	options = spot_sdp_default_options();
	options.verbose=1;

	sol=prog.minimize(-sum(slack),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		% Vertices_values=double(sol.eval(Vertices_values));
		w1=double(sol.eval(diff(V1,x)))
		w2=double(sol.eval(diff(V2,x)))
		w3=double(sol.eval(diff(V3,x)))
		w4=double(sol.eval(diff(V4,x)))
		delta_rho=double(sol.eval(delta_rho))
		V1=(sol.eval(V1));
		V2=(sol.eval(V2));
		V3=(sol.eval(V3));
		V4=(sol.eval(V4));
		new_v_values=double(sol.eval(this_vert_values))
		% rho=double(sol.eval(rho));
		% w=[w1,w2,w3,w4]./rho;
	else
		sol_OK=false;
	end
end


function [rho,Vertices_values,w,sol_OK]=diamond(x,xdot,df,method,varargin)
	% disp(method);
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	rho=varargin{1};
	[prog,Vertices_values] = prog.newPos(4);
	prog=prog.withPos(Vertices_values()-1e-6*ones(4,1));
	Lmonom = monomials(x,0:2);
	[prog,L] = prog.newSOSPoly(Lmonom,12);
	% the normals, without scaling (so that w is at the same scale as vertice values)
	w1=[-Vertices_values(2);Vertices_values(1)];
	w2=[-Vertices_values(2);-Vertices_values(3)];
	w3=[Vertices_values(4);-Vertices_values(3)];
	w4=[Vertices_values(4);Vertices_values(1)];
	% slack variables, pushing the solution into the interior of the feasible set
	[prog,slack]=prog.newPos(4);

	constraint1=[x(1);-x(2);-x(1)+x(2)-rho;];
	constraint2=[x(1);x(2);-x(1)-x(2)-rho;];
	constraint3=[x(2);-x(1);x(1)-x(2)-rho;];
	constraint4=[-x(1);-x(2);x(1)+x(2)-rho;];



	% polytope 1
	% V1dot=w1'*xdot;
	% constraint1=x(1);
	% constraint2=-x(2);
	% constraint3=-x(1)+x(2)-rho;
	prog=prog.withSOS(-slack(1)-V1dot+L(1:3)*constraint1);

	% polytope 2
	% V2dot=w2'*xdot;
	% constraint4=x(1);
	% constraint5=x(2);
	% constraint6=-x(1)-x(2)-rho;
	prog=prog.withSOS(-slack(2)-V2dot+L(4:6)*constraint2);

	% polytope 3
	% V3dot=w3'*xdot;
	% constraint7=x(2);
	% constraint8=-x(1);
	% constraint9=x(1)-x(2)-rho;
	prog=prog.withSOS(-slack(3)-V3dot+L(7:9)*constraint3);

	% polytope 4
	% V4dot=w4'*xdot;
	% constraint10=-x(1);
	% constraint11=-x(2);
	% constraint12=x(1)+x(2)-rho;
	prog=prog.withSOS(-slack(4)-V4dot+L(10:12)*constraint4);
	options = spot_sdp_default_options();
	options.verbose=1;

	sol=prog.minimize(-sum(slack),@spot_mosek,options);

	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		Vertices_values=double(sol.eval(Vertices_values));
		w1=double((sol.eval(w1)));
		w2=double(sol.eval(w2));
		w3=double(sol.eval(w3));
		w4=double(sol.eval(w4));
		rho=double(sol.eval(rho));
		% getting the right scale of the original w
		w=[w1,w2,w3,w4]./rho;

	else 
		sol_OK=false;
		% falls back to the last valid rho
		rho=varargin{1}/1.2;
		Vertices_values=varargin{2};
		w=varargin{3};
		figure(2)
		syms x1 x2;
		V= piecewise(x1<=0&-x2<=0&-x1+x2-rho<=0,w(:,1)'*[x1;x2],x1<=0&x2<=0&-x1-x2-rho<=0,w(:,2)'*[x1;x2],x2<=0&-x1<=0&x1-x2-rho<=0,w(:,3)'*[x1;x2],-x1<=0&-x2<=0&x1+x2-rho<=0,w(:,4)'*[x1;x2]);
		fsurf(1e3*V); hold on
		fcontour(V);
		figure(3)
		sym_xdot =[x2; -x1-x2.*(x1.^2-1)];
		Vdot=piecewise(x1<=0&-x2<=0&-x1+x2-rho<=0,w(:,1)'*sym_xdot,x1<=0&x2<=0&-x1-x2-rho<=0,w(:,2)'*sym_xdot,x2<=0&-x1<=0&x1-x2-rho<=0,w(:,3)'*sym_xdot,-x1<=0&-x2<=0&x1+x2-rho<=0,w(:,4)'*sym_xdot);
		fsurf(1e3*Vdot); hold on
	end
end



