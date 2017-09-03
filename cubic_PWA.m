function cubic_PWA()
checkDependency('spotless');
checkDependency('mosek');
x=msspoly('x',1);
xdot = [-x+x.^3];

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog,slack] = prog.newPos(1);

Vmonom = monomials(x,0:1);
[prog,V] = prog.newFreePoly(Vmonom);

Lmonom = monomials(x,0:2);
[prog,L1] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L1);

old_rho=1;
vert1=subs(V,x,[0]);
vert2=subs(V,x,[old_rho]);

prog=prog.withEqs(vert1);
prog=prog.withPos(vert2-1);

Vdot=diff(V,x)*xdot;
constraint1=-x(1)*(old_rho-x(1));
prog=prog.withSOS(-Vdot+L1*constraint1);

options = spot_sdp_default_options();
options.verbose=1;
sol=prog.minimize(0,@spot_mosek,options);
V=sol.eval(V)
% slack=sol.eval(slack)
L1=sol.eval(L1)
subs(L1,x,5)
end