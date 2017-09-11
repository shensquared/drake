function plots_stuff(x,xdot,V,all_V,last_rho,delta_rho)
        sum_rho=last_rho+delta_rho;
        constraint1=[x(1);-x(2);-x(1)+x(2)-sum_rho;-(-x(1)+x(2)-last_rho)];
        constraint2=[x(1);x(2);-x(1)-x(2)-sum_rho;-(-x(1)-x(2)-last_rho)];
        constraint3=[x(2);-x(1);x(1)-x(2)-sum_rho;-(x(1)-x(2)-last_rho)];
        constraint4=[-x(1);-x(2);x(1)+x(2)-sum_rho;-(x(1)+x(2)-last_rho)];
        figure(1)
        subplot(1,4,1)
        combined=[-x(1)+x(2)-sum_rho;-x(1)-x(2)-sum_rho;x(1)-x(2)-sum_rho;x(1)+x(2)-sum_rho;];
        [a,b]=meshgrid(-sum_rho:sum_rho/100:sum_rho,-sum_rho:sum_rho/100:sum_rho);
        regional=[min(dmsubs(combined,x,[a(:)';b(:)'])<=0)];
        z=(dmsubs(V,x,[a(:)';b(:)']));
        z=regional.*max(z,[],1);
        z(z==0) = NaN;
        z=reshape(z,size(a));
        surf(a,b,z);
        title('$$max(V)$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        hold on;
        subplot(1,4,2)
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
        subplot(1,4,3)
        surf(a,b,z);
        title('piecewise affine $$V$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        subplot(1,4,4)
        contour(a,b,z);
        title('piecewise affine $$V$$ contour' ,'interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        figure(2)
        % Vdot plot
        % subplot(1,5,5)
        w=diff(V,x);
        regional=[min(dmsubs(constraint1,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint2,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint3,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint4,x,[a(:)';b(:)'])<=0)];
        z=(regional).*(dmsubs(w*xdot,x,[a(:)';b(:)']));
        z=sum(z,1);
        z=reshape(z,size(a));
        z(z==0) = NaN;
        surf(a,b,z);
        title('piecewise $$\dot{V}$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
end