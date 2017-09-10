function [V,rho,all_V,sol_OK]=VanDerPol_PWA()
	do_plots=false;
	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);

	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	% xdot=[-x(1);-x(2)];
	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];

	rho=3e-4;
	delta_rho=3e-4;
	Vertices_values=zeros(4,1);
	w=zeros(2,4);
	sol_OK=true;
	level=0;
	[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,0,level,rho);
	level=1;
	while(sol_OK)
		disp(level);
		[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,rho,1,delta_rho,V,all_V);
		level=level+1;
    end

	% diamond_ring(x,xdot,rho,[7.6180;7.5818;7.6180;7.5818],.2);
% 	[rho,Vertices_values,w,sol_OK]=diamond_ring(x,xdot,rho,Vertices_values,1,do_plots);
% 	[rho,Vertices_values,w,sol_OK]=diamond_ring(x,xdot,rho,Vertices_values,.4,do_plots);

end


function [V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,varargin)
	last_rho=varargin{1};
	level=varargin{2};
	delta_rho=varargin{3};
	sum_rho=last_rho+delta_rho;
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);

	Vmonom = monomials(x,0:1);
	[prog,V] = prog.newFreePoly(Vmonom,4);
	
	Lmonom = monomials(x,0:2);
	[prog,L] = prog.newSOSPoly(Lmonom,16);

	vert0=[0;0];
	vert1=[0;last_rho];
	vert2=[-last_rho;0];
	vert3=[0;-last_rho];
	vert4=[last_rho;0];
	vert5=[0;last_rho+delta_rho];
	vert6=[-last_rho-delta_rho;0];
	vert7=[0;-last_rho-delta_rho];
	vert8=[last_rho+delta_rho;0];
	inner_verts=[vert1,vert2,vert3,vert4];
	outter_verts=[vert5,vert6,vert7,vert8];

	if level==0
		% all V zero at zero
		prog=prog.withEqs(subs(V,x,vert0)-zeros(4,1));
		% all V strictly positive at verts
		prog=prog.withPos((subs(V(1),x,vert5))-1e-6);
		prog=prog.withPos((subs(V(4),x,vert5))-1e-6);
		prog=prog.withPos((subs(V(1),x,vert6))-1e-6);
		prog=prog.withPos((subs(V(2),x,vert6))-1e-6);
		prog=prog.withPos((subs(V(2),x,vert7))-1e-6);
		prog=prog.withPos((subs(V(3),x,vert7))-1e-6);
		prog=prog.withPos((subs(V(3),x,vert8))-1e-6);
		prog=prog.withPos((subs(V(4),x,vert8))-1e-6);
	else
		last_V=varargin{4};
		all_V=varargin{5};
		last_max=max((dmsubs(last_V,x,inner_verts)),[],2);
		% evaluated at inner verts, last_V>=V
		prog=prog.withPos(last_max(1)-(subs(V(1),x,vert1)));
		prog=prog.withPos(last_max(1)-(subs(V(4),x,vert1)));
		prog=prog.withPos(last_max(2)-(subs(V(1),x,vert2)));
		prog=prog.withPos(last_max(2)-(subs(V(2),x,vert2)));
		prog=prog.withPos(last_max(3)-(subs(V(2),x,vert3)));
		prog=prog.withPos(last_max(3)-(subs(V(3),x,vert3)));
		prog=prog.withPos(last_max(4)-(subs(V(3),x,vert4)));
		prog=prog.withPos(last_max(4)-(subs(V(4),x,vert4)));
		% V at outter verts are bigger than inner verts
		prog=prog.withPos((subs(V(1),x,vert5-vert1)));
		prog=prog.withPos((subs(V(4),x,vert5-vert1)));
		prog=prog.withPos((subs(V(1),x,vert6-vert2)));
		prog=prog.withPos((subs(V(2),x,vert6-vert2)));
		prog=prog.withPos((subs(V(2),x,vert7-vert3)));
		prog=prog.withPos((subs(V(3),x,vert7-vert3)));
		prog=prog.withPos((subs(V(3),x,vert8-vert4)));
		prog=prog.withPos((subs(V(4),x,vert8-vert4)));
	end

	w=diff(V,x);
	% slack variables, pushing the solution into the interior of the feasible set
	[prog,slack]=prog.newPos(4);

	constraint1=[x(1);-x(2);-x(1)+x(2)-sum_rho;-(-x(1)+x(2)-last_rho)];
	constraint2=[x(1);x(2);-x(1)-x(2)-sum_rho;-(-x(1)-x(2)-last_rho)];
	constraint3=[x(2);-x(1);x(1)-x(2)-sum_rho;-(x(1)-x(2)-last_rho)];
	constraint4=[-x(1);-x(2);x(1)+x(2)-sum_rho;-(x(1)+x(2)-last_rho)];

	% polytope 1
	V1dot=w(1,:)*xdot*delta_rho;
	prog=prog.withSOS(-slack(1)-V1dot+L(1:4)'*constraint1);

	% polytope 2
	V2dot=w(2,:)*xdot*delta_rho;
	prog=prog.withSOS(-slack(2)-V2dot+L(5:8)'*constraint2);

	% polytope 3
	V3dot=w(3,:)*xdot*delta_rho;
	prog=prog.withSOS(-slack(3)-V3dot+L(9:12)'*constraint3);

	% polytope 4
	V4dot=w(4,:)*xdot*delta_rho;
	prog=prog.withSOS(-slack(4)-V4dot+L(13:16)'*constraint4);
	
	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(slack),@spot_mosek,options);

	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
        V=sol.eval(V);
        w=sol.eval(w);
        rho=sum_rho;
        figure(1)
        % subplot(1,5,1)
        clf
        combined=[-x(1)+x(2)-sum_rho;-x(1)-x(2)-sum_rho;x(1)-x(2)-sum_rho;x(1)+x(2)-sum_rho;];
        [a,b]=meshgrid(-sum_rho:sum_rho/50:sum_rho,-sum_rho:sum_rho/50:sum_rho);
        regional=[min(dmsubs(combined,x,[a(:)';b(:)'])<=0)];
        if level==0
        	all_V=V;
        	z=regional.*(dmsubs([V],x,[a(:)';b(:)']));
        else
        	all_V=[all_V;V];
        	z=regional.*(dmsubs(all_V,x,[a(:)';b(:)']));
        end
        z=max(z,[],1);
        z(z==0) = NaN;
        z=reshape(z,size(a));
        surf(a,b,z);
        title('$$max(V)$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        % subplot(1,5,2)
        figure(2)
        clf
        contour(a,b,z);
        title('$$ max(V)$$ contour','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        % true piecewise
        % regional=[min(dmsubs(constraint1,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint2,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint3,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint4,x,[a(:)';b(:)'])<=0)];
        % z=(regional).*(dmsubs(V,x,[a(:)';b(:)']));
        % z=max(z,[],1);
        % z=reshape(z,size(a));
        % z(z==0) = NaN;
        % subplot(1,5,3)
        % surf(a,b,z);
        % title('piecewise affine $$V$$','interpreter','latex','fontsize',20) 
        % xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        % ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        % subplot(1,5,4)
        % contour(a,b,z);
        % title('piecewise affine $$V$$ contour' ,'interpreter','latex','fontsize',20) 
        % xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        % ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        % % Vdot plot
        % subplot(1,5,5)
        % V=sol.eval(V);
        % w=sol.eval(w);
        % regional=[min(dmsubs(constraint1,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint2,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint3,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint4,x,[a(:)';b(:)'])<=0)];
        % z=(regional).*(dmsubs(w*xdot,x,[a(:)';b(:)']));
        % z=sum(z,1);
        % z=reshape(z,size(a));
        % z(z==0) = NaN;
        % surf(a,b,z);
        % title('piecewise $$\dot{V}$$','interpreter','latex','fontsize',20) 
        % xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        % ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	else 
		sol_OK=false;
		% falls back to the last valid rho
		rho=varargin{1}/1.2;
		Vertices_values=varargin{2};
		w=varargin{3};
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



function [V,rho]=discontinuous_diamond(x,xdot,do_plots,varargin)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	rho=varargin{1};
	level=varargin{2};
	last_V=varargin{3};

	Vmonom = monomials(x,0:1);
	[prog,V] = prog.newFreePoly(Vmonom,4);

	Lmonom = monomials(x,0:2);
	[prog,L] = prog.newSOSPoly(Lmonom,12);

	vert0=[0;0];
	vert1=[0;rho];
	vert2=[-rho;0];
	vert3=[-rho;0];
	vert4=[0;-rho];
	vert5=[0;-rho];
	vert6=[rho;0];
	vert7=[rho;0];
	vert8=[0;rho];



	prog=prog.withEqs(subs(V,x,vert0)-zeros(4,1));
% all V strictly positive at verts
	prog=prog.withPos((subs(V(1),x,vert1))-1e-6);
	prog=prog.withPos((subs(V(1),x,vert2))-1e-6);
	prog=prog.withPos((subs(V(2),x,vert3))-1e-6);
	prog=prog.withPos((subs(V(2),x,vert4))-1e-6);
	prog=prog.withPos((subs(V(3),x,vert5))-1e-6);
	prog=prog.withPos((subs(V(3),x,vert6))-1e-6);
	prog=prog.withPos((subs(V(4),x,vert7))-1e-6);
	prog=prog.withPos((subs(V(4),x,vert8))-1e-6);


% the normals, without scaling (so that w is at the same scale as vertice values)
	w=diff(V,x);
	% w1=[-Vertices_values(2);Vertices_values(1)];
	% w2=[-Vertices_values(2);-Vertices_values(3)];
	% w3=[Vertices_values(4);-Vertices_values(3)];
	% w4=[Vertices_values(4);Vertices_values(1)];
	% slack variables, pushing the solution into the interior of the feasible set
	
	[prog,slack]=prog.newPos(4);

	constraint1=[x(1);-x(2);-x(1)+x(2)-rho;];
	constraint2=[x(1);x(2);-x(1)-x(2)-rho;];
	constraint3=[x(2);-x(1);x(1)-x(2)-rho;];
	constraint4=[-x(1);-x(2);x(1)+x(2)-rho;];

	% polytope 1
	V1dot=w(1,:)*xdot*rho;
	prog=prog.withSOS(-slack(1)-V1dot+L(1:3)'*constraint1);

	% polytope 2
	V2dot=w(2,:)*xdot*rho;
	prog=prog.withSOS(-slack(2)-V2dot+L(4:6)'*constraint2);

	% polytope 3
	V3dot=w(3,:)*xdot*rho;
	prog=prog.withSOS(-slack(3)-V3dot+L(7:9)'*constraint3);

	% polytope 4
	V4dot=w(4,:)*xdot*rho;
	prog=prog.withSOS(-slack(4)-V4dot+L(10:12)'*constraint4);
	
	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(slack),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
        V=sol.eval(V);
        w=sol.eval(w);
        figure(1)
        subplot(1,5,1)
        combined=[-x(1)+x(2)-rho;-x(1)-x(2)-rho;x(1)-x(2)-rho;x(1)+x(2)-rho;];
        [a,b]=meshgrid(-rho:rho/100:rho,-rho:rho/100:rho);
        regional=[min(dmsubs(combined,x,[a(:)';b(:)'])<=0)];
        z=regional.*(dmsubs(V,x,[a(:)';b(:)']));
        z=max(z,[],1);
        z(z==0) = NaN;
        z=reshape(z,size(a));
        surf(a,b,z);
        title('$$max(V)$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        hold on;
        subplot(1,5,2)
        contour(a,b,z);
        title('$$ max(V)$$ contour','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        hold on;
        % true piecewise
        regional=[min(dmsubs(constraint1,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint2,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint3,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint4,x,[a(:)';b(:)'])<=0)];
        z=(regional).*(dmsubs(V,x,[a(:)';b(:)']));
        z=max(z,[],1);
        z=reshape(z,size(a));
        z(z==0) = NaN;
        subplot(1,5,3)
        surf(a,b,z);
        title('piecewise affine $$V$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        subplot(1,5,4)
        contour(a,b,z);
        title('piecewise affine $$V$$ contour' ,'interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        % Vdot plot
        subplot(1,5,5)
        V=sol.eval(V);
        w=sol.eval(w);
        regional=[min(dmsubs(constraint1,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint2,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint3,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint4,x,[a(:)';b(:)'])<=0)];
        z=(regional).*(dmsubs(w*xdot,x,[a(:)';b(:)']));
        z=sum(z,1);
        z=reshape(z,size(a));
        z(z==0) = NaN;
        surf(a,b,z);
        title('piecewise $$\dot{V}$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	else 
		sol_OK=false;
		% falls back to the last valid rho
		rho=varargin{1}/1.2;
		Vertices_values=varargin{2};
		w=varargin{3};
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



function [rho,Vertices_values,w,sol_OK]=continuous_diamond_ring(x,xdot,df,last_rho,last_vertice_values,delta_rho,do_plots)
	prog = spotsosprog;
	prog = prog.withIndeterminate(x);
	[prog,this_vert_values] = prog.newPos(4);
	prog=prog.withPos(this_vert_values-1e-6*last_rho*ones(4,1));
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

	prog=prog.withPos(-(subs(V1,x,vert1)-last_vertice_values(1)));
	prog=prog.withPos(-(subs(V1,x,vert2)-last_vertice_values(2)));
	prog=prog.withPos(-(subs(V2,x,vert2)-last_vertice_values(2)));
	prog=prog.withPos(-(subs(V2,x,vert3)-last_vertice_values(3)));
	prog=prog.withPos(-(subs(V3,x,vert3)-last_vertice_values(3)));
	prog=prog.withPos(-(subs(V3,x,vert4)-last_vertice_values(4)));
	prog=prog.withPos(-(subs(V4,x,vert4)-last_vertice_values(4)));
	prog=prog.withPos(-(subs(V4,x,vert1)-last_vertice_values(1)));


	prog=prog.withPos((subs(V1,x,vert5)-last_vertice_values(1)));
	prog=prog.withPos((subs(V1,x,vert6)-last_vertice_values(2)));
	prog=prog.withPos((subs(V2,x,vert6)-last_vertice_values(2)));
	prog=prog.withPos((subs(V2,x,vert7)-last_vertice_values(3)));
	prog=prog.withPos((subs(V3,x,vert7)-last_vertice_values(3)));
	prog=prog.withPos((subs(V3,x,vert8)-last_vertice_values(4)));
	prog=prog.withPos((subs(V4,x,vert8)-last_vertice_values(4)));
	prog=prog.withPos((subs(V4,x,vert5)-last_vertice_values(1)));


	% prog=prog.withEqs(subs(V1,x,vert5)-this_vert_values(1));
	% prog=prog.withEqs(subs(V1,x,vert5)-this_vert_values(1));
	% prog=prog.withEqs(subs(V1,x,vert6)-this_vert_values(2));
	% prog=prog.withEqs(subs(V2,x,vert6)-this_vert_values(2));
	% prog=prog.withEqs(subs(V2,x,vert7)-this_vert_values(3));
	% prog=prog.withEqs(subs(V3,x,vert7)-this_vert_values(3));
	% prog=prog.withEqs(subs(V3,x,vert8)-this_vert_values(4));
	% prog=prog.withEqs(subs(V4,x,vert8)-this_vert_values(4));
	% prog=prog.withEqs(subs(V4,x,vert5)-this_vert_values(1));
	
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
	V1dot=diff(V1,x)*xdot*sum_rho;
	V2dot=diff(V2,x)*xdot*sum_rho;
	V3dot=diff(V3,x)*xdot*sum_rho;
	V4dot=diff(V4,x)*xdot*sum_rho;
	prog=prog.withSOS(-slack(1)-V1dot+L(1:4)'*constraint1);
	prog=prog.withSOS(-slack(2)-V2dot+L(5:8)'*constraint2);
	prog=prog.withSOS(-slack(3)-V3dot+L(9:12)'*constraint3);
	prog=prog.withSOS(-slack(4)-V4dot+L(13:16)'*constraint4);

	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(0),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		% Vertices_values=double(sol.eval(Vertices_values));
		w1=double(sol.eval(diff(V1,x)));
		w2=double(sol.eval(diff(V2,x)));
		w3=double(sol.eval(diff(V3,x)));
		w4=double(sol.eval(diff(V4,x)));
		w=[w1,w2,w3,w4];
		delta_rho=double(sol.eval(delta_rho));
		b1=double(sol.eval(subs(V1,x,zeros(2,1))));
		b2=double(sol.eval(subs(V2,x,zeros(2,1))));
		b3=double(sol.eval(subs(V3,x,zeros(2,1))));
		b4=double(sol.eval(subs(V4,x,zeros(2,1))));
		% b1=0;
		% b2=0;
		% b3=0;
		% b4=0;

		V1=(sol.eval(V1));
		V2=(sol.eval(V2));
		V3=(sol.eval(V3));
		V4=(sol.eval(V4));

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





% 