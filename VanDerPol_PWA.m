function VanDerPol_PWA()
	flags=struct();
	flags.method='diamond';
	% flags.method='square';
	flags.scaling=false;
	flags.do_plots=true;
	flags.verbose=true;
	flags.debug=true;
	flags.quiver=false;

	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);

	xdot=[-x(1);-x(2)];
% 	xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	% xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];
	
	if flags.quiver
		plot_quiver(x,xdot,last_rho)
    end
    
    last_rho=0;
	delta_rho=3e-2;
	last_V=zeros(4,1);
	all_V=[];
	sol_OK=true;
	[last_V,last_rho,all_V,sol_OK]=flat_square_rings(x,xdot,last_rho,delta_rho,last_V,all_V,flags);
end
