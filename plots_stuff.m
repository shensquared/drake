function plots_stuff(x,xdot,V,all_V,last_rho,delta_rho,flags)
        sum_rho=last_rho+delta_rho;

        switch flags.method
        case 'diamond'        
                constraint1=[x(1);-x(2);-x(1)+x(2)-sum_rho;-(-x(1)+x(2)-last_rho)];
                constraint2=[x(1);x(2);-x(1)-x(2)-sum_rho;-(-x(1)-x(2)-last_rho)];
                constraint3=[x(2);-x(1);x(1)-x(2)-sum_rho;-(x(1)-x(2)-last_rho)];
                constraint4=[-x(1);-x(2);x(1)+x(2)-sum_rho;-(x(1)+x(2)-last_rho)];
                combined=[-x(1)+x(2)-sum_rho;-x(1)-x(2)-sum_rho;x(1)-x(2)-sum_rho;x(1)+x(2)-sum_rho;];
        case 'square'
                constraint1=[-sum_rho-x(1);x(1)+x(2);x(1)-x(2);-(-last_rho-x(1))];
                constraint2=[-sum_rho-x(2);x(2)-x(1);x(1)+x(2);-(-last_rho-x(2))];
                constraint3=[x(2)-x(1);-x(1)-x(2);x(1)-sum_rho;-(x(1)-last_rho)];
                constraint4=[x(2)-sum_rho;-x(1)-x(2);x(1)-x(2);-(x(2)-last_rho)];
                combined=[-sum_rho-x(1);-sum_rho-x(2);x(1)-sum_rho;x(2)-sum_rho];
        end
        [a,b]=meshgrid(-sum_rho:sum_rho/100:sum_rho,-sum_rho:sum_rho/100:sum_rho);

        figure(1)
        clf
        % true piecewise
        regional=[min(dmsubs(constraint1,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint2,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint3,x,[a(:)';b(:)'])<=0);min(dmsubs(constraint4,x,[a(:)';b(:)'])<=0)];
        z=(regional).*(dmsubs(V,x,[a(:)';b(:)']));
        z=sum(z,1);
        z=reshape(z,size(a));
        z(z==0) = NaN;
        subplot(1,2,1)
        surf(a,b,z);
        title('piecewise affine $$V$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        subplot(1,2,2)
        contour(a,b,z,10);
        title('piecewise affine $$V$$ contour' ,'interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        figure(2)
        % Vdot plot
        subplot(1,2,1)
        w=diff(V,x);
        regional=[min(dmsubs(constraint1,x,[a(:)';b(:)'])<0);min(dmsubs(constraint2,x,[a(:)';b(:)'])<0);min(dmsubs(constraint3,x,[a(:)';b(:)'])<0);min(dmsubs(constraint4,x,[a(:)';b(:)'])<0)];
        z=(regional).*(dmsubs(w*xdot,x,[a(:)';b(:)']));
        z=sum(z,1);
        z=reshape(z,size(a));
        z(z==0) = NaN;
        surf(a,b,z);
        title('piecewise $$\dot{V}$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
        subplot(1,2,2)
        contour(a,b,z,10);
        title('piecewise affine $$\dot{V}$$ contour' ,'interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)

end