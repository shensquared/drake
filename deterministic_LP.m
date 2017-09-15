function deterministic_LP()
	sol_OK=false;

	rho=.9;
	resolution=18;
	row_verts=2*resolution+1;
	num_tris=(2*resolution)^2;
	num_verts=(row_verts)^2;


	[a,b]=meshgrid(-rho:rho/resolution:rho,-rho:rho/resolution:rho);
	[centroid_a,centroid_b]=meshgrid((-3+1/resolution)*rho/3:rho/resolution:(3-2/resolution)*rho/3,(-3+1/resolution)*rho/3:rho/resolution:(3-2/resolution)*rho/3);
	[centroid_c,centroid_d]=meshgrid((-3+2/resolution)*rho/3:rho/resolution:(3-1/resolution)*rho/3,(-3+2/resolution)*rho/3:rho/resolution:(3-1/resolution)*rho/3);

	centroid_a=reshape(centroid_a,[num_tris,1]);
	centroid_b=reshape(centroid_b,[num_tris,1]);
	centroid_c=reshape(centroid_c,[num_tris,1]);
	centroid_d=reshape(centroid_d,[num_tris,1]);
	% clf
	% scatter(reshape(a,[num_verts,1]),reshape(b,[num_verts,1])); hold on
	% scatter(centroid_a,centroid_b); hold on
	% scatter(centroid_c,centroid_d); hold on

	x=msspoly('x',2);
	% xdot = [-2*x(1); -2*x(2)];
	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];

	df=diff(xdot,x);
	A=[];
	A1=dmsubs(df(:,1),x,[centroid_a(:)';centroid_b(:)']);
	A2=dmsubs(df(:,2),x,[centroid_a(:)';centroid_b(:)']);
	for i =1:num_tris
		A=[A,A1(:,i),A2(:,i)];
	end
	B=[];
	B1=dmsubs(df(:,1),x,[centroid_c(:)';centroid_d(:)']);
	B2=dmsubs(df(:,2),x,[centroid_c(:)';centroid_d(:)']);
	for i =1:num_tris
		B=[B,B1(:,i),B2(:,i)];
	end
	vert_values=find_V(x,xdot,A,B,rho,resolution,a,b);
end

function vert_values=find_V(x,xdot,A,B,rho,resolution,a,b)
	row_verts=2*resolution+1;
	num_tris=(2*resolution)^2;
	num_verts=(row_verts)^2;
	prog=spotsosprog;
	prog=prog.withIndeterminate(x);
	[prog,vert_values]=prog.newPos(num_verts);
	% [prog,slacks]=prog.newPos(num_tris);
	% slacks=reshape(slacks,[2*resolution,2*resolution]);

	% [prog,slacks_B]=prog.newPos(num_tris);
	% slacks_B=reshape(slacks_B,[2*resolution,2*resolution]);

% vertices values constrains
	vert_values=reshape(vert_values,[row_verts,row_verts]);
	% vert_value=0 at the origin
	prog=prog.withEqs(vert_values(resolution+1,resolution+1));
	% the exterior values are the same, all equal to 1
	prog=prog.withEqs(vert_values(1,:)-1);
	prog=prog.withEqs(vert_values(row_verts,:)-1);
	prog=prog.withEqs(vert_values(:,1)-1);
	prog=prog.withEqs(vert_values(:,row_verts)-1);

	column_diff=vert_values(:,2)-vert_values(:,1);
	row_diff=vert_values(2,:)-vert_values(1,:);

	for i=2:row_verts-1
		column_diff=[column_diff,vert_values(:,i+1)-vert_values(:,i)];
		row_diff=[row_diff;vert_values(i+1,:)-vert_values(i,:)];
	end

% % lower triangles
	w1=column_diff(2:row_verts,:);
	w2=row_diff(:,1:row_verts-1);
	w1=reshape(w1,num_tris,1);
	w2=reshape(w2,num_tris,1);
	xrep=repmat(x,[num_tris,1]);
	Vdot=[w1,w2]*A*xrep;
	Vdot=reshape(Vdot,[2*resolution,2*resolution]);
	for i=1:2*resolution
		for j=1:2*resolution
			prog=prog.withPos(-(subs(Vdot(i,j),x,[a(i,j);b(i,j)])));
			prog=prog.withPos(-(subs(Vdot(i,j),x,[a(i+1,j);b(i+1,j)])));
			prog=prog.withPos(-(subs(Vdot(i,j),x,[a(i+1,j+1);b(i+1,j+1)])));
		end
	end
% upper triangless
	n1=column_diff(1:row_verts-1,:);
	n2=row_diff(:,2:row_verts);
	n1=reshape(n1,num_tris,1);
	n2=reshape(n2,num_tris,1);

	Vdot_B=[n1,n2]*B*xrep;
	Vdot_B=reshape(Vdot_B,[2*resolution,2*resolution]);
	% for i=1:2*resolution
	% 	for j=1:2*resolution
	% 		% verts=[(a(i,j);b(i,j));;(a(i+1,j);b(i+1,j))];
	% 		prog=prog.withPos(-(subs(Vdot_B(i,j),x,[a(i,j);b(i,j)])));
	% 		prog=prog.withPos(-(subs(Vdot_B(i,j),x,[a(i,j+1);b(i,j+1)])));
	% 		prog=prog.withPos(-(subs(Vdot_B(i,j),x,[a(i+1,j+1);b(i+1,j+1)])));

	% 	end
	% end

	% slacks=reshape(slacks,[num_tris,1]);
	% slacks_B=reshape(slacks_B,[num_tris,1]);

	
	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(0,@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		% sol.eval(slacks)
		values=double(sol.eval(vert_values))
		tri = delaunay(a,b);
		trisurf(tri,a,b,values);
		title('piecewise $$\dot{V}$$','interpreter','latex','fontsize',20) 
		xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
		ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	end

end



