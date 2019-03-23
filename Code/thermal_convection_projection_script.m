[U,V]=thermal_convection_projection();
%%
figure(1);
t = 0.5;
subplot(2,2,1);contourf(U'); colorbar;
x = linspace(0,1,size(U,1));y = linspace(0,1,size(U,2)); [x,y]=meshgrid(x,y);
dx = 1/size(U,1);
u = (y.^2-y)*t;
m = 20;
errorU = min(abs((U'-u)./(u+eps)), m*dx);
errorU = abs(U'-u);
subplot(2,2,2); contourf(errorU); colorbar;
subplot(2,2,3);contourf(V'); colorbar;
x = linspace(0,1,size(V,1));y = linspace(0,1,size(V,2)); [x,y]=meshgrid(x,y);
v = (x.^2-x)*t;
errorV = min(abs((V'-v)./(v+eps)),m*dx);
errorV = abs(V'-v);
subplot(2,2,4); contourf(errorV); colorbar;
disp([norm(errorU,2)/numel(errorU),norm(errorV,2)/numel(errorV)]);
