function [V,rho,all_V,sol_OK]=VanDerPol_PWA()
	flags=struct();
	flags.method='diamond';
% 	flags.method='discontinuous';
	% flags.method='square'
	flags.do_plots=true;
	flags.verbose=true;
	flags.debug=true;
	flags.quiver=false;

	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);

	% xdot=[-x(1);-x(2)];
	xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
% 	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];

	last_rho=.8;
    
	if flags.quiver
		plot_quiver(x,xdot,last_rho)
    end
    
	delta_rho=.1;
	last_V=zeros(4,1);
	all_V=[];
	sol_OK=true;

	switch flags.method
	case 'diamond'
		[last_V,last_rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,last_rho,delta_rho,last_V,all_V,flags)
	case 'square'
		[last_V,last_rho,all_V,sol_OK]=flat_square_rings(x,xdot,last_rho,delta_rho,last_V,all_V,flags)
	case 'discontinuous'
		[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,last_rho,delta_rho,last_V,all_V,flags)
	end

	% while(sol_OK)
	% 	[last_V,last_rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,level,last_rho,delta_rho,last_V,all_V,do_plots,verbose)
	% 	level=level+1;
	% end

	% [V,rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,level,rho,delta_rho,V,all_V,do_plots)
	% [rho,Vertices_values,w,sol_OK]=cont_diamond_ring(x,xdot,rho,zeros(4,1),delta_rho,do_plots)
end
