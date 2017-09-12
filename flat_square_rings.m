function [V,rho,all_V,sol_OK]=flat_square_rings(x,xdot,last_rho,delta_rho,last_V,all_V,flags)
	do_plots=flags.do_plots;
	verbose=flags.verbose;
	debug_flag=flags.debug;
	method=flags.method;

	sum_rho=last_rho+delta_rho;
	sum_rho_order=10^(floor(log10(sum_rho)));	
	delta_rho_order=10^(floor(log10(sum_rho)));	

	prog = spotsosprog;
	prog = prog.withIndeterminate(x);

	% the normals
	[prog,w]=prog.newPos(1);
	% prog=prog.withPos(w-1e-3);
	% slack variables, pushing the solution into the interior of the feasible set
	[prog,slack]=prog.newPos(4);
	
	Lmonom = monomials(x,1:2);
	LwConstMonom=monomials(x,1:2);

	[prog,this_flat_value]=prog.newPos(1);
	prog=prog.withPos(this_flat_value-1e-6*delta_rho);

	switch flags.method
	case 'diamond'
		vert1=[0;last_rho];
		vert2=[-last_rho;0];
		vert3=[0;-last_rho];
		vert4=[last_rho;0];
		vert5=[0;last_rho+delta_rho];
		vert6=[-last_rho-delta_rho;0];
		vert7=[0;-last_rho-delta_rho];
		vert8=[last_rho+delta_rho;0];
		V1=-w*x(1)+w*x(2);
		V2=-w*x(1)-w*x(2);
		V3=w*x(1)-w*x(2);
		V4=w*x(1)+w*x(2);

	case 'square'
		vert1=[-last_rho;last_rho];
		vert2=[-last_rho;-last_rho];
		vert3=[last_rho;-last_rho];
		vert4=[last_rho;last_rho];
		vert5=[-sum_rho;sum_rho];
		vert6=[-sum_rho;-sum_rho];
		vert7=[sum_rho;-sum_rho];
		vert8=[sum_rho;sum_rho];
		V1=-w*x(1);
		V2=-w*x(2);
		V3=w*x(1);
		V4=w*x(2);
	end
	inner_verts=[vert1,vert2,vert3,vert4];
	outter_verts=[vert5,vert6,vert7,vert8];

	V1dot=diff(V1,x)*xdot;
	V2dot=diff(V2,x)*xdot;
	V3dot=diff(V3,x)*xdot;
	V4dot=diff(V4,x)*xdot;
	
	if last_rho==0
		[prog,L] = prog.newSOSPoly(Lmonom,8);
		[prog,LwConst] = prog.newSOSPoly(LwConstMonom,4);
		[prog,l]=prog.newPos(12);
		switch flags.method
		case 'diamond'
			constraint1=[x(1);-x(2);-x(1)+x(2)-sum_rho;];
			constraint2=[x(1);x(2);-x(1)-x(2)-sum_rho;];
			constraint3=[x(2);-x(1);x(1)-x(2)-sum_rho];
			constraint4=[-x(1);-x(2);x(1)+x(2)-sum_rho];
		case 'square'
			constraint1=[-sum_rho-x(1);x(1)+x(2);x(1)-x(2)];
			constraint2=[-sum_rho-x(2);x(2)-x(1);x(1)+x(2)];
			constraint3=[x(2)-x(1);-x(1)-x(2);x(1)-sum_rho];
			constraint4=[x(2)-sum_rho;-x(1)-x(2);x(1)-x(2)];
		end

		if flags.scaling
			constraint1=diag([1,1,1/sum_rho_order])*constraint1;
			constraint2=diag([1,1,1/sum_rho_order])*constraint2;
			constraint3=diag([1,1,1/sum_rho_order])*constraint3;
			constraint4=diag([1,1,1/sum_rho_order])*constraint4;
		end

		prog=prog.withSOS((-slack(1)-(sum_rho_order)^2*V1dot+[L(1:2)',LwConst(1)']*constraint1+l(1:3)'*constraint1));
		prog=prog.withSOS((-slack(2)-(sum_rho_order)^2*V2dot+[L(3:4)',LwConst(2)']*constraint2+l(4:6)'*constraint2));
		prog=prog.withSOS((-slack(3)-(sum_rho_order)^2*V3dot+[L(5:6)',LwConst(3)']*constraint3+l(7:9)'*constraint3));
		prog=prog.withSOS((-slack(4)-(sum_rho_order)^2*V4dot+[L(7:8)',LwConst(4)']*constraint4+l(10:12)'*constraint4));	

	else
		[prog,L] = prog.newSOSPoly(Lmonom,8);
		[prog,LwConst] = prog.newSOSPoly(LwConstMonom,8);
		[prog,l]=prog.newPos(16);
		
		switch flags.method
		case 'diamond'
			constraint1=[x(1);-x(2);-x(1)+x(2)-sum_rho;-(-x(1)+x(2)-last_rho)];
			constraint2=[x(1);x(2);-x(1)-x(2)-sum_rho;-(-x(1)-x(2)-last_rho)];
			constraint3=[x(2);-x(1);x(1)-x(2)-sum_rho;-(x(1)-x(2)-last_rho)];
			constraint4=[-x(1);-x(2);x(1)+x(2)-sum_rho;-(x(1)+x(2)-last_rho)];
		case 'square'
			constraint1=[-sum_rho-x(1);x(1)+x(2);x(1)-x(2);-(-last_rho-x(1))];
	        constraint2=[-sum_rho-x(2);x(2)-x(1);x(1)+x(2);-(-last_rho-x(2))];
	        constraint3=[x(2)-x(1);-x(1)-x(2);x(1)-sum_rho;-(x(1)-last_rho)];
	        constraint4=[x(2)-sum_rho;-x(1)-x(2);x(1)-x(2);-(x(2)-last_rho)];
		end

		if flags.scaling
			constraint1=diag([1,1,1/sum_rho_order,1/delta_rho_order])*constraint1;
			constraint2=diag([1,1,1/sum_rho_order,1/delta_rho_order])*constraint2;
			constraint3=diag([1,1,1/sum_rho_order,1/delta_rho_order])*constraint3;
			constraint4=diag([1,1,1/sum_rho_order,1/delta_rho_order])*constraint4;
		end
		prog=prog.withSOS((-slack(1)-(sum_rho_order)^2*V1dot+[L(1:2)',LwConst(1:2)']*constraint1)+l(1:4)'*constraint1);
		prog=prog.withSOS((-slack(2)-(sum_rho_order)^2*V2dot+[L(3:4)',LwConst(3:4)']*constraint2)+l(5:8)'*constraint2);
		prog=prog.withSOS((-slack(3)-(sum_rho_order)^2*V3dot+[L(5:6)',LwConst(5:6)']*constraint3)+l(9:12)'*constraint3);
		prog=prog.withSOS((-slack(4)-(sum_rho_order)^2*V4dot+[L(7:8)',LwConst(7:8)']*constraint4)+l(13:16)'*constraint4);	

	end

	options = spot_sdp_default_options();
	options.verbose=verbose;
	sol=prog.minimize(sum(-slack),@spot_mosek,options);

	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
        V=[sol.eval(V1);sol.eval(V2);sol.eval(V3);sol.eval(V4)];
        this_flat_value=sol.eval(this_flat_value);
        if last_rho==0
        	all_V=V;
        else
        	all_V=[all_V;V];
        end
        plots_stuff(x,xdot,V,all_V,last_rho,delta_rho,flags);        
        rho=sum_rho;
        if debug_flag
        	disp('L')
        	L=sol.eval(L)
        	sol.eval(LwConst)
        	sol.eval(l)
        	% disp('V')
        	% sol.eval(V) 
        	disp('Vdot')
        	sol.eval(diff(V,x)*xdot)
        	disp('flat_value')
        	this_flat_value
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

