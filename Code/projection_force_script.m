[U,V,f1,f2]=projection_force(425,425);
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
%% computing errors
i = 1;
K = floor(30*1.7.^[0:4]);
errors_u=[]; errors_u=[];
for k = K
    [U_k,V_k]=projection_force(k,k);
    x = ceil(linspace(1,size(U,1),size(U_k,1)));
    y = ceil(linspace(1,size(U,2),size(U_k,2)));
    errors_u(i) = norm((U(x,y)-U_k),2)/sqrt(numel(U_k));
    x = ceil(linspace(1,size(V,1),size(V_k,1)));
    y = ceil(linspace(1,size(V,2),size(V_k,2)));
    errors_v(i) = norm((V(x,y)-V_k),2)/sqrt(numel(V_k));
    i = i + 1;
end
% Estimating order
eu_ratio = errors_u(2:end)./errors_u(1:end-1);
ev_ratio = errors_v(2:end)./errors_v(1:end-1);
kratio = K(2:end)./K(1:end-1);
log(eu_ratio)./log(kratio)