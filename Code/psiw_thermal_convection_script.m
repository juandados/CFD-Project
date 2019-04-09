%% Glycerin Setup
[xx,yy,Us,Vs,Ts,Cs,Ws,dt,tend,tS,Re,Pr,Ra,Pe] = psiw_thermal_convection(20,'g',48,8);
%%
figure(2);
u = Us(:,:,end);v=Vs(:,:,end);T=Ts(:,:,end);C=Cs(:,:,end);w=Ws(:,:,end);
%u = Us(:,:,1);v=Vs(:,:,1);T=Ts(:,:,1);C=Cs(:,:,1);w=Ws(:,:,1);
Len = sqrt(u.^2+v.^2)+eps;
dx = xx(1,2)-xx(1,1);
dy = yy(2,1)-yy(1,1);
lx = xx(1,end); ly = yy(end,1);
pcolor(0:dx:lx,0:dy:ly,T); shading interp;
contourf(0:dx:lx,0:dy:ly,T);
colorbar; hold on;
%quiver(xx,yy,(u./Len),(v./Len),0.6,'k-');
contour(xx,yy,w,7,'w-');
axis equal, axis([0 lx 0 ly]);
hold off, axis equal, axis([0 lx 0 ly]);
%title(sprintf('Numerical Results Glycerine Benchmark',Re,tS));
drawnow
print('Aspect Ratio 6:1','-dpdf', '-fillpage');
%% Plotting resuts Air;
figure(1);
u = Us(:,:,end);v=Vs(:,:,end);T=Ts(:,:,end);C=Cs(:,:,end);w=Ws(:,:,end);
u = Us(:,:,1);v=Vs(:,:,1);T=Ts(:,:,1);C=Cs(:,:,1);w=Ws(:,:,1);
Len = sqrt(u.^2+v.^2)+eps;
dx = xx(1,2)-xx(1,1);
dy = yy(2,1)-yy(1,1);
lx = xx(1,end); ly = yy(end,1);
contourf(0:dx:lx,0:dy:ly,T);
%pcolor(0:dx:lx,0:dy:ly,T); shading interp; 
colorbar; hold on;
quiver(xx,yy,(u./Len),(v./Len),0.6,'k-');
contour(xx,yy,w,20,'w-');
axis equal, axis([0 lx 0 ly]);
hold off, axis equal, axis([0 lx 0 ly]);
title(sprintf('Re = %0.1g   t = %0.2g',Re,tS));
drawnow
%% simulate for the errors
clear all;
powers=[3:-1:3];
UU = struct(); VV = struct(); TT = struct(); CC = struct(); WW = struct();
XX = struct(); YY = struct();
for k = powers
    [xx,yy,Us,Vs,Ts,Cs,Ws,dt,tend,tS,Re,Pr,Ra,Pe] = psiw_thermal_convection(30,'a',4*2^k,1*2^k);
    key = strcat('k',num2str(k));
    UU = setfield(UU,key,Us);
    VV = setfield(VV,key,Vs);
    TT = setfield(TT,key,Ts);
    CC = setfield(CC,key,Cs);
    WW = setfield(WW,key,Ws);
    XX = setfield(XX,key,xx);
    YY = setfield(YY,key,yy);
end
save(strrep(strcat(sprintf('air 4-1 Re: %0.2f Pr: %0.2f Ra %0.2f Pe %0.2f #',Re,Pr,Ra,Pe),datestr(datetime())),'.','_'))
%% compute errors
% load manually the .mat file with the data if variables not avaliable
close all;
powers=[2:5];
kf = powers(end);
U = getfield(UU,strcat('k',num2str(kf)));
V = getfield(VV,strcat('k',num2str(kf)));
T = getfield(TT,strcat('k',num2str(kf)));
C = getfield(CC,strcat('k',num2str(kf)));
W = getfield(WW,strcat('k',num2str(kf)));
Uerrors = []; Verrors = []; Terrors = [];
for k = powers(1:end-1)
    key = strcat('k',num2str(k));
    Uk = getfield(UU,key);
    Vk = getfield(VV,key);
    Tk = getfield(TT,key);
    Ck = getfield(CC,key);
    Wk = getfield(WW,key);
    events_count = size(Uk,3);
    for event = 1:events_count
        %figure();
        relative_error = min(abs(U(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event)-Uk(1:end-1,1:end-1,event))./(eps+abs(U(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event))),10);
        %contourf(relative_error);colorbar; title(['U k',num2str(k),'event',num2str(event)]);
        Uerrors(k,event) = norm(relative_error,2)/sqrt(numel(relative_error));
        relative_error = min(abs(V(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event)-Vk(1:end-1,1:end-1,event))./(eps+abs(V(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event))),10);
        Verrors(k,event) = norm(relative_error,2)/sqrt(numel(relative_error));
        relative_error = min(abs(T(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event)-Tk(1:end-1,1:end-1,event))./(eps+abs(T(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event))),10);
        Terrors(k,event) = norm(relative_error,2)/sqrt(numel(relative_error));
        relative_error = min(abs(C(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event)-Ck(1:end-1,1:end-1,event))./(eps+abs(C(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event))),10);
        Cerrors(k,event) = norm(relative_error,2)/sqrt(numel(relative_error));
        relative_error = min(abs(W(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event)-Wk(1:end-1,1:end-1,event))./(eps+abs(W(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event))),10);
        Werrors(k,event) = norm(relative_error,2)/sqrt(numel(relative_error));
    end
end
%% Plotting solution at steady state
key = strcat('k',num2str(5));
Uk = getfield(UU,key);
figure(91); contourf(Uk(:,:,end)); colorbar;
% key = strcat('k',num2str(3));
% Uk = getfield(UU,key);
% figure(94); contourf(Uk(:,:,end-4)); colorbar;

%% Convergence plots
figure()
mesh(linspace(2,4,size(Werrors,1)-1)',linspace(0,5e3,size(Werrors,2)-1)',Werrors(2:end,2:end)');
zlabel('Relative error %');ylabel('time');xlabel('-log(\Delta x)');
%%
print('Relative error V%','-dpdf');
%% finding Ps
errors = Cerrors(3:end,end);
logs=log2(errors(1:end-1)./errors(2:end));
p = mean(logs(end-2:end))
%% plotting results
% load manually the .mat file with the data if variables not avaliable
close all;
for k = powers(1:end)
    key = strcat('k',num2str(k));
    Uk = getfield(UU,key);
    Vk = getfield(VV,key);
    Tk = getfield(TT,key);
    Ck = getfield(CC,key);
    Wk = getfield(WW,key);
    XXk = getfield(XX,key);
    YYk = getfield(YY,key);
    events_count = size(Uk,3);
    figure();
    for event = 1:events_count
        subplot(events_count,1,event)
        u = Uk(:,:,event);v=Vk(:,:,event);T=Tk(:,:,event);
        C=Ck(:,:,event);w=Wk(:,:,event);xx=XXk; yy=YYk;
        Len = sqrt(u.^2+v.^2)+eps;
        dx = xx(1,2)-xx(1,1);
        dy = yy(2,1)-yy(1,1);
        lx = xx(1,end); ly = yy(end,1);
        contourf(0:dx:lx,0:dy:ly,T);
        %pcolor(0:dx:lx,0:dy:ly,T); shading interp; 
        colorbar; hold on;
        quiver(xx,yy,(u./Len),(v./Len),0.6,'k-');
        contour(xx,yy,w,5,'w-');
        axis equal, axis([0 lx 0 ly]);
        hold off, axis equal, axis([0 lx 0 ly]);
        title(sprintf('Re = %0.1g   t = %0.2g',Re,event*tS/events_count));
    end
    print('planktonVT','-dpdf', '-fillpage');
    figure();
    for event = 1:events_count
        subplot(events_count,1,event);
        u = Uk(:,:,event);v=Vk(:,:,event);T=Tk(:,:,event);
        C=Ck(:,:,event);w=Wk(:,:,event);xx=XXk; yy=YYk;
        Len = sqrt(u.^2+v.^2)+eps;
        dx = xx(1,2)-xx(1,1);
        dy = yy(2,1)-yy(1,1);
        lx = xx(1,end); ly = yy(end,1);
        contourf(0:dx:lx,0:dy:ly,C);
        %pcolor(0:dx:lx,0:dy:ly,T); shading interp; 
        colorbar; hold on;
        quiver(xx,yy,(u./Len),(v./Len),0.6,'k-');
        contour(xx,yy,w,5,'w-');
        axis equal, axis([0 lx 0 ly]);
        hold off, axis equal, axis([0 lx 0 ly]);
        title(sprintf('C Re = %0.1g   t = %0.2g',Re,event*tS/events_count));
    end
    print('planktonC','-dpdf', '-fillpage');
end 
%% Graphic for poster
% load manually the .mat file with the data if variables not avaliable
close all;
key = strcat('k','6');
Uk = getfield(UU,key);
Vk = getfield(VV,key);
Tk = getfield(TT,key);
Ck = getfield(CC,key);
Wk = getfield(WW,key);
XXk = getfield(XX,key);
YYk = getfield(YY,key);
events = [3,4,12];
for m = 1:3 
    fig=figure();
    event = events(m);
    u = Uk(:,:,event);v=Vk(:,:,event);T=Tk(:,:,event);
    C=Ck(:,:,event);w=Wk(:,:,event);xx=XXk; yy=YYk;
    ss = 1:5:size(xx,1);
    pp = 1:5:size(xx,2);
    Len = sqrt(u(ss,pp).^2+v(ss,pp).^2)+eps;
    dx = xx(1,2)-xx(1,1);
    dy = yy(2,1)-yy(1,1);
    lx = xx(1,end); ly = yy(end,1);
    contourf(0:dx:lx,0:dy:ly,T);
    %pcolor(0:dx:lx,0:dy:ly,T); shading interp; 
    %colorbar; 
    hold on;
    %quiver(xx(ss,pp),yy(ss,pp),(u(ss./Len),(v./Len),0.6,'k-');
    quiver(xx(ss,pp),yy(ss,pp),u(ss,pp)./Len,v(ss,pp)./Len,0.6,'k-');
    axis equal, axis([0 lx 0 ly]);
    hold off, axis equal, axis([0 lx 0 ly]);
    title(sprintf('Velocity and Temperature for t = %0.2f s', event*tS/events_count));
    %subplot(2,3,m+3);
    orient(fig,'landscape');
    print(['T',num2str(m)],'-dpdf', '-fillpage');
    fig=figure();
    contourf(0:dx:lx,0:dy:ly,C);
    %pcolor(0:dx:lx,0:dy:ly,T); shading interp; 
    %colorbar; 
    hold on;
    contour(xx,yy,w,5,'w-','LineWidth',1);
    axis equal, axis([0 lx 0 ly]);
    hold off, axis equal, axis([0 lx 0 ly]);
    title(['Concetration and Stream Function for t = ', num2str(event*tS/events_count),' s']);
    orient(fig,'landscape');
    print(['C',num2str(m)],'-dpdf', '-fillpage');
end
close all;

%% Get cells dependency
% from file: air Re: 4365_30 Pr: 0_72 Ra 300142_85 Pe 0_07 #31-Mar-2019 16:00:46
% load manually the .mat file with the data if variables not avaliable
close all;
kk = 0;
for k = [4,6]
    key = strcat('k',num2str(k));
    Uk = getfield(UU,key);
    Vk = getfield(VV,key);
    Tk = getfield(TT,key);
    Ck = getfield(CC,key);
    Wk = getfield(WW,key);
    XXk = getfield(XX,key);
    YYk = getfield(YY,key);
    events_count = size(Uk,3);
    event = 5;
    kk = kk + 1;
    subplot(2,1,kk);
    u = Uk(:,:,event);v=Vk(:,:,event);T=Tk(:,:,event);
    C=Ck(:,:,event);w=Wk(:,:,event);xx=XXk; yy=YYk;
    Len = sqrt(u.^2+v.^2)+eps;
    dx = xx(1,2)-xx(1,1);
    dy = yy(2,1)-yy(1,1);
    lx = xx(1,end); ly = yy(end,1);
    contourf(0:dx:lx,0:dy:ly,T);
    %pcolor(0:dx:lx,0:dy:ly,T); shading interp; 
    colorbar;hold on;
    xs=[1:(size(xx,2)-1)/64:size(xx,2)];
    ys=[1:(size(xx,1)-1)/16:size(xx,1)];
    quiver(xx(ys,xs),yy(ys,xs),(u(ys,xs)./Len(ys,xs)),(v(ys,xs)./Len(ys,xs)),0.6,'k-');
    contour(xx,yy,w,5,'w-','LineWidth',1);
    axis equal, axis([0 lx 0 ly]);
    hold off, axis equal, axis([0 lx 0 ly]);
    % title('t=1000s')
    ue
    colorbar;
end 
print('cells','-dpdf', '-fillpage');

%% Plankton case
% from file: load('plackton Re: 184269_66 Pr: 7_56 Ra 96726349664_99 Pe 69653_93 #01-Apr-2019 00:12:22.mat')
% run the print section an then save
print('plankton','-dpdf', '-fillpage');
%% Get Dependency on time step:
%% run air setup for nx=2^3 dt=1e-1 for:
% k = k*0.1, k=10*k keeping others as nominal for air
% Re = Re*0.1, Re=10*k
% Pr = Pr*0.1, Pr=10*k
[xx,yy,Us,Vs,Ts,Cs,Ws,dt,tend,tS,Re,Pr,Ra,Pe] = psiw_thermal_convection(20,'a',48,8);
figure();
key = strcat('k','3');
u = Us(:,:,end);v=Vs(:,:,end);T=Ts(:,:,end);
w=Ws(:,:,end);
Len = sqrt(u.^2+v.^2)+eps;
dx = xx(1,2)-xx(1,1);
dy = yy(2,1)-yy(1,1);
lx = xx(1,end); ly = yy(end,1);
colorbar; hold on;
contourf(0:dx:lx,0:dy:ly,T);
quiver(xx,yy,(u./Len),(v./Len),0.6,'k-');
%contour(xx,yy,w,5,'w-');
axis equal, axis([0 lx 0 ly]);
hold off, axis equal, axis([0 lx 0 ly]);
print('param dependence air re x0_1','-dpdf', '-fillpage');
%% Get Dependency on time step:
% run air setup for nx=2^3 dt=1e-1
% run air setup for nx=2^3 dt=1e-2
figure();
key = strcat('k','3');
Uk = getfield(UU,key);
Vk = getfield(VV,key);
Tk = getfield(TT,key);
Ck = getfield(CC,key);
Wk = getfield(WW,key);
XXk = getfield(XX,key);
YYk = getfield(YY,key);
u = Uk(:,:,end);v=Vk(:,:,end);T=Tk(:,:,end);
xx=XXk; yy=YYk; w=Wk(:,:,end);
Len = sqrt(u.^2+v.^2)+eps;
dx = xx(1,2)-xx(1,1);
dy = yy(2,1)-yy(1,1);
lx = xx(1,end); ly = yy(end,1);
colorbar; hold on;
contourf(0:dx:lx,0:dy:ly,T);
quiver(xx,yy,(u./Len),(v./Len),0.6,'k-');
contour(xx,yy,w,5,'w-');
axis equal, axis([0 lx 0 ly]);
hold off, axis equal, axis([0 lx 0 ly]);
print('time step structure dependence dt 1e-1','-dpdf', '-fillpage');