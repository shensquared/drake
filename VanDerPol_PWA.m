function [V,rho,all_V,sol_OK]=VanDerPol_PWA()
	flags=struct();

	flags.method='flat';
% 	flags.method='discontinuous';
	flags.do_plots=false;
	flags.verbose=true;
	flags.debug=true;
	flags.quiver=false;

	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);

	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];

	% xdot=[-x(1);-x(2)];
	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];

	last_rho=0;
	delta_rho=2e-2;
	last_V=zeros(4,1);
	all_V=[];
	sol_OK=true;

	switch flags.method
	case 'flat'
		[last_V,last_rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,last_rho,delta_rho,last_V,all_V,flags)
	case 'discontinuous'
		[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,last_rho,delta_rho,last_V,all_V,flags)
	end
	if flags.quiver
		plot_quiver(x,xdot,last_rho)
	% while(sol_OK)
	% 	[last_V,last_rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,level,last_rho,delta_rho,last_V,all_V,do_plots,verbose)
	% 	level=level+1;
	% end

	% [V,rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,level,rho,delta_rho,V,all_V,do_plots)
	% [rho,Vertices_values,w,sol_OK]=cont_diamond_ring(x,xdot,rho,zeros(4,1),delta_rho,do_plots)
end
