function ran_tri(x,xdot,flags)
	containment=false;
	% generating samples containning the origin
	samples=zeros(3,2);
	while (~containment)
		samples=randn(3,2);
		containment=checkContainment(x,samples);
	end
	% disp(samples)
	V=oneLevelV(x,samples)

	V1dot=diff(V(1),x)*xdot;
	V2dot=diff(V(2),x)*xdot;
	V3dot=diff(V(3),x)*xdot;
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);

	Lmonom = monomials(x,0:2);
	[prog,L] = prog.newSOSPoly(Lmonom,9);

	constraint1=[V(1)-1;V(2)-V(1);V(3)-V(1)];
	constraint2=[V(2)-1;V(1)-V(2);V(3)-V(2)];
	constraint3=[V(3)-1;V(1)-V(3);V(2)-V(3)];

	[prog,slack]=prog.newPos(3);
	% [prog,scalings]=prog.newPos(3);

	prog=prog.withSOS(-slack(1)-V1dot+[L(1:3)']*constraint1);
	prog=prog.withSOS(-slack(2)-V2dot+[L(4:6)']*constraint2);
	prog=prog.withSOS(-slack(3)-V3dot+[L(7:9)']*constraint3);

	options = spot_sdp_default_options();
	options.verbose=0;
	sol=prog.minimize(-sum(slack),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		disp('L')
		sol.eval(L)
		tri_plots(x,xdot,samples,V,flags)
	else
		error('derivative condition can not be satisified')
	end

end

function containment=checkContainment(x,samples)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	[prog,lambda]=prog.newPos(3);
	prog=prog.withEqs(sum(lambda)-1);
	prog=prog.withEqs(lambda'*samples);
	options = spot_sdp_default_options();
	options.verbose=0;
	sol=prog.minimize(0,@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		containment=true;
	else
		containment=false;
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