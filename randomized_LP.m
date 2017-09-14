function [samples,V,rho,sol_OK]=randomized_LP(x,xdot,flags,rho)
	containment=false;
	sol_OK=false;
	% generating two linearly indepdenet samples
	linear_depedent=false;
	samples=zeros(2,2);
	while (~linear_depedent)
		samples=randn(2,2);
		linear_depedent=rank(samples==2);
	end
	centroid=2/3*[sum(samples,1);sum(samples,2)];

	% disp(samples)

	while(~sol_OK)
		[rho,V,sol_OK]=line_search_rho(x,xdot,rho,V,samples,flags);
		rho=rho/1.2
	end


end

function [V,level_values]=find_V(x,xdot,samples,centroid,flags)
	prog=spotsosprog;
	prog=prog.withIndeterminate(x);
	Vmonom=monomials(x,1:1);
	[prog,V]=prog.newFreePoly(Vmonom);

	% zero at zero
	prog=prog.withEqs(subs(V,x,zeros(2,1)));
	% strictly positive at the other vertices
	prog=prog.withPos(subs(V,x,samples));

	df=diff(xdot,x);%nbyn
	A=subs(df,x,centroid);

	partialV=diff(V,x);  %1byn
	% vdot<0 at the vertices, with f approxiamted at the centroid
	prog=prog.withPos(-(subs(partialV*A*x,x,samples)));
	
	options = spot_sdp_default_options();
	options.verbose=0;
	sol=prog.minimize(-sum(0),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		sol.eval(L)
		tri_plots(x,xdot,samples,V,flags)
	else
		sol_OK=false;
		if rho<1e-8
			error('derivative condition can not be satisified')
		end
	end	
end




function [rho,V,sol_OK]=line_search_rho(x,xdot,rho,V,samples,flags)
	V1dot=diff(V(1),x)*xdot;
	V2dot=diff(V(2),x)*xdot;
	V3dot=diff(V(3),x)*xdot;
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);

	Lmonom = monomials(x,0:2);
	[prog,L] = prog.newSOSPoly(Lmonom,9);

	constraint1=[V(1)-rho;V(2)-V(1);V(3)-V(1)];
	constraint2=[V(2)-rho;V(1)-V(2);V(3)-V(2)];
	constraint3=[V(3)-rho;V(1)-V(3);V(2)-V(3)];

	[prog,slack]=prog.newPos(3);
	% [prog,scalings]=prog.newPos(3);

	prog=prog.withSOS(-slack(1)-V1dot+[L(1:3)']*constraint1);
	prog=prog.withSOS(-slack(2)-V2dot+[L(4:6)']*constraint2);
	prog=prog.withSOS(-slack(3)-V3dot+[L(7:9)']*constraint3);

	options = spot_sdp_default_options();
	options.verbose=0;
	sol=prog.minimize(-sum(slack),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		disp('L')
		sol.eval(L)
		tri_plots(x,xdot,samples,V,flags)
	else
		sol_OK=false;
		if rho<1e-8
			error('derivative condition can not be satisified')
		end
	end	
end


function V=oneLevelV(x,samples)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	Vmonom = monomials(x,1:1);
	[prog,V] = prog.newFreePoly(Vmonom,3);
	prog=prog.withEqs(subs(V(1),x,samples(1,:)')-1);
	prog=prog.withEqs(subs(V(2),x,samples(1,:)')-1);
	prog=prog.withEqs(subs(V(2),x,samples(2,:)')-1);
	prog=prog.withEqs(subs(V(3),x,samples(2,:)')-1);
	prog=prog.withEqs(subs(V(3),x,samples(3,:)')-1);
	prog=prog.withEqs(subs(V(1),x,samples(3,:)')-1);

	options = spot_sdp_default_options();
	options.verbose=0;
	sol=prog.minimize(0,@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		V=sol.eval(V);
	else
		error('can not find a valid tuple of candidate V')
    end
end