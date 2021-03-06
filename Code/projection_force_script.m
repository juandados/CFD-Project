[U,V,f1,f2]=projection_force(2^9,2^9);
%%
figure(1);
t = 0.5e-1;
subplot(2,2,1);contourf(U'); colorbar;
x = linspace(0,1,size(U,1));y = linspace(0,1,size(U,2)); [x,y]=meshgrid(x,y);
dx = 1/size(U,1);
u = (y.^2-y)*t;
m = 200;
errorU = min(abs(U'-u)./(abs(u)+eps), 1e-1);
%errorU = abs(U'-u);
subplot(2,2,2); contourf(errorU); colorbar;
subplot(2,2,3);contourf(V'); colorbar;
x = linspace(0,1,size(V,1));y = linspace(0,1,size(V,2)); [x,y]=meshgrid(x,y);
v = (x.^2-x)*t;
errorV = min(abs(V'-v)./(abs(v)+eps),1e-1);
%errorV = abs(V'-v);
subplot(2,2,4); contourf(errorV); colorbar;
disp([norm(errorU,2)/numel(errorU),norm(errorV,2)/numel(errorV)]);
%% computing errors
i = 1;
K = floor(2.^[4:8]);
errors_u=[]; errors_v=[];
for k = K
    [U_k,V_k]=projection_force(k,k);
    % computing U accuracy
    intercell_distance_y = size(U,2)/size(U_k,2);
    intercell_distance_x = (size(U,1)+1)/(size(U_k,1)+1);
    x = [intercell_distance_x:intercell_distance_x:size(U,1)];
    y = [intercell_distance_y/2:intercell_distance_y:size(U,2)-intercell_distance_y/2];
    UU = (U(x,y)+U(x,y+1))/2;
    errors_u(i) = norm((UU-U_k),2)/sqrt(numel(U_k));
    % computing V accuracy
    intercell_distance_y = (size(V,2)+1)/(size(V_k,2)+1);
    intercell_distance_x = size(V,1)/size(V_k,1);
    x = [intercell_distance_x/2:intercell_distance_x:size(V,1)-intercell_distance_x/2];
    y = [intercell_distance_y:intercell_distance_y:size(V,2)];
    VV = (V(x,y)+V(x+1,y))/2;
    errors_v(i) = norm((VV-V_k),2)/sqrt(numel(V_k));
    i = i + 1;
end
% Estimating order
eu_ratio = errors_u(2:end)./errors_u(1:end-1);
ev_ratio = errors_v(2:end)./errors_v(1:end-1);
kratio = K(2:end)./K(1:end-1);
disp(['convergence rate space for U: ',num2str(-log(eu_ratio)./log(kratio))])
disp(['convergence rate space for V: ',num2str(-log(ev_ratio)./log(kratio))])