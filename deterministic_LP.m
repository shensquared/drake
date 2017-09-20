function deterministic_LP()
	findV=true;
	pause_sec=0;

	rho=.9;
	resolution=10;

	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);

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
	
	% xdot = [-2*x(1); -2*x(2)];
	% xdot = [-x(1)+x(1)^3; -x(2)+x(2)^3];
	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];
	A=[];Ax=[];bias_A=[];B=[];Bx=[];bias_B=[];
	for i =1:num_tris
		shifted_xdot=subs(xdot,x,[x(1)+centroid_a(i);x(2)+centroid_b(i)]);
		current_df=diff(shifted_xdot,x);
		current_centroid=[centroid_a(i)';centroid_b(i)'];
		current_bias=subs(xdot,x,current_centroid);
		bias_A=[bias_A,current_bias];
		A1=dmsubs(current_df(:,1),x,zeros(2,1));
		A2=dmsubs(current_df(:,2),x,zeros(2,1));
		current_A=[A1,A2];
		A=[A,A1,A2];
		Ax=[Ax,subs(current_A,x,zeros(2,1))];

		% shifted_xdot=subs(xdot,x,[x(1)-centroid_c(i);x(2)-centroid_d(i)]);
		% current_df=diff(shifted_xdot,x);
		% current_centroid=[centroid_c(i)';centroid_d(i)'];
		% current_bias=subs(xdot,x,current_centroid);
		% bias_B=[bias_B,current_bias];
		% B1=dmsubs(current_df(:,1),x,zeros(2,1));
		% B2=dmsubs(current_df(:,2),x,zeros(2,1));
		% current_B=[B1,B2];
		% B=[B,B1,B2];
		% Bx=[Bx,subs(current_B,x,zeros(2,1))];
	end

	% figure('units','normalized','outerposition',[0 0 1 1])
	% clf
	% scatter(reshape(a,[num_verts,1]),reshape(b,[num_verts,1]));
	% 	xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
	% 	ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	% axis equal tight
	% hold on
	% pause(pause_sec);
	% scatter(centroid_a,centroid_b); hold on
	% 	xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
	% 	ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	% pause(pause_sec);
	% % quiver(centroid_a,centroid_b,Ax(1,:)',Ax(2,:)'); hold on
	% % 	xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
	% % 	ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	% % pause(pause_sec);
	% scatter(centroid_c,centroid_d); hold on
	% 	xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
	% 	ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	% pause(pause_sec);
	% % quiver(centroid_c,centroid_d,Bx(1,:)',Bx(2,:)'); hold on
	% % 	xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
	% % 	ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	% % quiver(a,b,-a+a.^3,-b+b.^3);
	% quiver(a,b,-b,-(-a-b.*(a.^2-1)));

	if findV
		vert_values=find_V(x,xdot,A,B,rho,resolution,a,b,centroid_a,centroid_b,centroid_c,centroid_d,bias_A,bias_B);
	end
end

function vert_values=find_V(x,xdot,A,B,rho,resolution,a,b,centroid_a,centroid_b,centroid_c,centroid_d,bias_A,bias_B)
	row_verts=2*resolution+1;
	num_tris=(2*resolution)^2;
	num_verts=(row_verts)^2;
	prog=spotsosprog;
	prog=prog.withIndeterminate(x);
	[prog,vert_values]=prog.newPos(num_verts);
	[prog,slacks]=prog.newPos(num_tris);
	slacks=reshape(slacks,[2*resolution,2*resolution]);
	% [prog,slacks_B]=prog.newPos(num_tris);
	% slacks_B=reshape(slacks_B,[2*resolution,2*resolution]);
	[prog,ext_value]=prog.newPos(1);
	prog=prog.withEqs(ext_value-10);
% vertices values constrains
	vert_values=reshape(vert_values,[row_verts,row_verts]);
	vert_values=flip(vert_values',1);
	% vert_value=0 at the origin
	prog=prog.withEqs(vert_values(resolution+1,resolution+1));
	% the exterior values are the same, all >ext_value
	prog=prog.withEqs(vert_values(1,:)-ext_value);
	prog=prog.withEqs(vert_values(row_verts,:)-ext_value);
	prog=prog.withEqs(vert_values(:,1)-ext_value);
	prog=prog.withEqs(vert_values(:,row_verts)-ext_value);
	prog=prog.withPos(ext_value*ones(size(vert_values(2:end,2:end)))-vert_values(2:end,2:end));
% partialVpartialx
	column_diff=flip(vert_values(:,2)-vert_values(:,1))';
	row_diff=(vert_values(end-1,:)-vert_values(end,:))';
	for i=2:row_verts-1
		column_diff=[column_diff;flip(vert_values(:,i+1)-vert_values(:,i))'];
		row_diff=[row_diff,(vert_values(end-i,:)-vert_values(end-i+1,:))'];
	end
	% column_diff=reshape(column_diff,[row_verts*(row_verts-1),1]);
% % % lower triangles
	w1=column_diff(:,1:row_verts-1);
	w2=row_diff(1:row_verts-1,:);
	w1=reshape(w1,num_tris,1);
	w2=reshape(w2,num_tris,1);
    w=[w1,w2];


	Vdot=w(1,:)*A(:,1:2)*x;
	bias_w=w(1,:)*bias_A(:,1);
	trueVdotA=w(1,:)*xdot;
	for i=2:num_tris
		Vdot=[Vdot;w(i,:)*A(:,2*i-1:2*i)*x];
		bias_w=[bias_w;w(i,:)*bias_A(:,i)];
		trueVdotA=[trueVdotA;w(i,:)*xdot];
	end
	Vdot=reshape(Vdot,[2*resolution,2*resolution]);
	bias_w=reshape(bias_w,[2*resolution,2*resolution]);
	trueVdotA=reshape(trueVdotA,[2*resolution,2*resolution]);

% upper triangless
	% n1=column_diff(2:row_verts,:);
	% n2=row_diff(:,2:row_verts);
	% n1=reshape(n1,num_tris,1);
	% n2=reshape(n2,num_tris,1);
 %    n=[n1,n2];
	% Vdot_B=n(1,:)*B(:,1:2)*x;
	% bias_w_B=n(1,:)*bias_B(:,1);
	% trueVdotB=n(1,:)*xdot;
	% for i=2:num_tris
	% 	Vdot_B=[Vdot_B;n(i,:)*B(:,2*i-1:2*i)*x];
	% 	bias_w_B=[bias_w_B;n(i,:)*bias_B(:,i)];
	% 	trueVdotB=[trueVdotB;n(i,:)*xdot];
 %    end
	% Vdot_B=reshape(Vdot_B,[2*resolution,2*resolution]);
	% bias_w_B=reshape(bias_w_B,[2*resolution,2*resolution]);
	% trueVdotB=reshape(trueVdotB,[2*resolution,2*resolution]);

	a=reshape(a,[num_verts,1]);
	b=reshape(b,[num_verts,1]);
	grids=[a';b'];
	grids=reshape(grids,[2,row_verts,row_verts]);

	centroid_a=reshape(centroid_a,[2*resolution,2*resolution]);
	centroid_b=reshape(centroid_b,[2*resolution,2*resolution]);
	centroid_c=reshape(centroid_c,[2*resolution,2*resolution]);
	centroid_d=reshape(centroid_d,[2*resolution,2*resolution]);

% Vdot<0 at grids
	for i=1:2*resolution
		for j=1:2*resolution
			% lower tri
			disp('lower')
			Vdotvert1=subs(Vdot(i,j),x,[grids(:,i,j)-[centroid_a(i,j);centroid_b(i,j)]])+bias_w(i,j)
			true1=subs(trueVdotA(i,j),x,[grids(:,i,j)])
			Vdotvert2=subs(Vdot(i,j),x,[grids(:,i,j+1)-[centroid_a(i,j);centroid_b(i,j)]])+bias_w(i,j)
			true2=subs(trueVdotA(i,j),x,[grids(:,i,j+1)])
			Vdotvert3=subs(Vdot(i,j),x,[grids(:,i+1,j)-[centroid_a(i,j);centroid_b(i,j)]])+bias_w(i,j)
			true3=subs(trueVdotA(i,j),x,[grids(:,i+1,j)])
			prog=prog.withPos(-slacks(i,j)-Vdotvert1);
			prog=prog.withPos(-slacks(i,j)-Vdotvert2);
			prog=prog.withPos(-slacks(i,j)-Vdotvert3);
			% upper tri
			% disp('upper')
			% Vdot_Bvert1=subs(Vdot_B(i,j),x,[grids(:,i+1,j)-[centroid_c(i,j);centroid_d(i,j)]])+bias_w_B(i,j)
			% trueB1=subs(trueVdotB(i,j),x,[grids(:,i+1,j)])
			% Vdot_Bvert2=subs(Vdot_B(i,j),x,[grids(:,i+1,j+1)-[centroid_c(i,j);centroid_d(i,j)]])+bias_w_B(i,j)
			% trueB2=subs(trueVdotB(i,j),x,[grids(:,i+1,j+1)])
			% Vdot_Bvert3=subs(Vdot_B(i,j),x,[grids(:,i,j+1)-[centroid_c(i,j);centroid_d(i,j)]])+bias_w_B(i,j)
			% trueB3=subs(trueVdotB(i,j),x,[grids(:,i,j+1)])
			% prog=prog.withPos(-slacks_B(i,j)-Vdot_Bvert1);
			% prog=prog.withPos(-slacks_B(i,j)-Vdot_Bvert2);
			% prog=prog.withPos(-slacks_B(i,j)-Vdot_Bvert3);
		end
	end


% solving stage
	slacks=reshape(slacks,[num_tris,1]);
	% slacks_B=reshape(slacks_B,[num_tris,1]);
	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(0,@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		% sol.eval(slacks)
		figure(25)
		values=double(sol.eval(vert_values))
		tri = delaunay(a,b);
		trisurf(tri,a,b,values);
		title('piecewise affine $${V}$$','interpreter','latex','fontsize',20) 
		xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
		ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	end
end



