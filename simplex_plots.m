% w=randn(3,2);
[a,b]=meshgrid(-5:.05:5,-5:.05:5);
V1=w(1,1)*a+w(1,2)*b;
V2=w(2,1)*a+w(2,2)*b;
V3=w(3,1)*a+w(3,2)*b;

% V1(V1>1)=NaN;
% V2(V2>1)=NaN;
% V3(V3>1)=NaN;

% V1(V1<0)=NaN;
% V2(V2<0)=NaN;
% V3(V3<0)=NaN;

% figure(1)
% clf
% surf(a,b,V1), hold on;
% surf(a,b,V2), hold on;
% surf(a,b,V3), hold on;

% V4=max(V1,V2,V3)
% surf(a,b,V4)
figure(2)
clf
% contour(a,b,V1,[0 0]), hold on;
% contour(a,b,V2,[0 0]), hold on;
% contour(a,b,V3,[0 0]), hold on;

[C1,h1]=contour(a,b,V1,[1 1]), hold on;
h1.LineColor='r';
[C2,h2]=contour(a,b,V2,[1 1]), hold on;
h2.LineColor='b';
[C3,h3]=contour(a,b,V3,[1 1]), hold on;
h3.LineColor='g'
% figure(3)
% clf
contour(a,b,V1-V2,[1 1]), hold on;
contour(a,b,V2-V3,[1 1]), hold on;
% contour(a,b,V3-V1,[1 1]), hold on;
contour(a,b,V1-V2,[-1 -1]), hold on;
contour(a,b,V2-V3,[-1 -1]), hold on;
% contour(a,b,V3-V1,[-1 -1]), hold on;
contour(a,b,V1-V2,[0 0]), hold on;
contour(a,b,V2-V3,[0 0]), hold on;
% contour(a,b,V3-V1,[0 0]), hold on;