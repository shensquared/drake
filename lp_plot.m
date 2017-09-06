
function lp_plot()
x1=-5:.01:5;
y1=-5:.01:5;
[X1,Y1]=meshgrid(x1,y1);
z1=abs(X1)+abs(Y1);

% surf(X1,Y1,z1);hold on
contour(X1,Y1,z1);

end
