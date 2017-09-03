function V=VanDerPol_PWA()
checkDependency('spotless');
checkDependency('mosek');
x=msspoly('x',2);
xdot = [x(2); -x(1)-x(2).*(x(1).^2-1)];
df = [zeros(2,1),[0 1; -1-2*x(1)*x(2), -(x(1)^2-1)]];

prog = spotsosprog;
prog = prog.withIndeterminate(x);
Vmonom = monomials(x,0:1);
[prog,V] = prog.newFreePoly(Vmonom);

[prog,slack] = prog.newPos(4);

old_rho=2e-4;
vert1=subs(V,x,[0;0]);
vert2=subs(V,x,[0;old_rho]);
vert3=subs(V,x,[old_rho;0]);
vert4=subs(V,x,[old_rho;old_rho]);
prog=prog.withEqs(vert1);
prog=prog.withPos(vert2);
prog=prog.withPos(vert3);
prog=prog.withPos(vert4);

Lmonom = monomials(x,0:2);
[prog,L1] = prog.newFreePoly(Lmonom);
[prog,L2] = prog.newFreePoly(Lmonom);
[prog,L3] = prog.newFreePoly(Lmonom);
[prog,L4] = prog.newFreePoly(Lmonom);

prog = prog.withSOS(L1);
prog = prog.withSOS(L2);
prog = prog.withSOS(L3);
prog = prog.withSOS(L4);

Vdot=diff(V,x)*xdot;
constraint1=-x(1);
constraint2=-x(2);
constraint3=-(old_rho-x(1));
constraint4=-(old_rho-x(2));

constraint5=-constraint1*constraint3;
constraint6=-constraint2*constraint4;
prog=prog.withSOS(-Vdot+L1*constraint1+L2*constraint2+L3*constraint3+L4*constraint4-slack(4));

% prog=prog.withSOS(-Vdot+L1*constraint5+L2*constraint6-slack(4));

options = spot_sdp_default_options();
options.verbose=1;
sol=prog.minimize(-vert4,@spot_mosek,options);
V=sol.eval(V);
end


function rho=find_boundary(V,x,options,old_rho)
prog = spotsosprog;
prog = prog.withIndeterminate(x);
Vmonom = monomials(x,0:1);
[prog,V] = prog.newFreePoly(Vmonom);

Lmonom = monomials(x,0:2);
[prog,L1] = prog.newFreePoly(Lmonom);
[prog,L2] = prog.newFreePoly(Lmonom);

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