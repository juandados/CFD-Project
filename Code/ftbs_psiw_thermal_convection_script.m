[xx,yy,Us,Vs,Ts,Cs,Ws,dt,tend,tS,Re,Pr,Ra,Pe] = ftbs_psiw_thermal_convection(10,'a',65,17);
%%
[u,v,T,w,ii,jj,dx,dy,b,g,xx,yy,Us,Vs,Ts,Cs,Ws,dt,tend,tS,Re,Pr,Ra,Pe]= ftbs_psiw_thermal_convection(10,'a',65,17);
%%  
r(jj,ii) = w(jj,ii) - dt * ...
      ((u(jj,ii)./(abs(u(jj,ii))+eps)).*u(jj,ii) .* (3*w(jj,ii)-4*w(jj,ii-1) + w(jj,ii-2)) / (2*dx) ...
	+ (v(jj,ii)./(abs(v(jj,ii))+eps)).*v(jj,ii) .* (3*w(jj,ii)-4*w(jj-1,ii) + w(jj-2,ii)) / (2*dy) ) ...
	+ (1/Re)*dt*(w(jj,ii+1) - 2*w(jj,ii) + w(jj,ii-1)) / dx^2 ...
	+ (1/Re)*dt*(w(jj+1,ii) - 2*w(jj,ii) + w(jj-1,ii)) / dy^2 ...
    + dt*b*g*(3*T(jj,ii)-4*T(jj,ii-1) + T(jj,ii-2)) / (2*dx);
%% Plotting resuts Air;
figure(1);
u = Us(:,:,end);v=Vs(:,:,end);T=Ts(:,:,end);C=Cs(:,:,end);w=Ws(:,:,end);
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
%% computing the errors
% simulate
clear all;
powers = [3:7];
UU = struct(); VV = struct(); TT = struct(); CC = struct(); WW = struct();
XX = struct(); YY = struct();
for k = powers
    [xx,yy,Us,Vs,Ts,Cs,Ws,dt,tend,tS,Re,Pr,Ra,Pe] = ftbs_psiw_thermal_convection(3,'a',4*2^k,1*2^k);
    key = strcat('k',num2str(k));
    UU = setfield(UU,key,Us);
    VV = setfield(VV,key,Vs);
    TT = setfield(TT,key,Ts);
    CC = setfield(CC,key,Cs);
    WW = setfield(WW,key,Ws);
    XX = setfield(XX,key,xx);
    YY = setfield(YY,key,yy);
end
save(strrep(strcat(sprintf('air Re: %0.2f Pr: %0.2f Ra %0.2f Pe %0.2f #',Re,Pr,Ra,Pe),datestr(datetime())),'.','_'))
%% compute errors
% load manually the .mat file with the data if variables not avaliable
kf = powers(end);
U = getfield(UU,strcat('k',num2str(kf)));
V = getfield(VV,strcat('k',num2str(kf)));
T = getfield(TT,strcat('k',num2str(kf)));
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
        relative_error = (U(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event)-Uk(1:end-1,1:end-1,event))/abs(eps+U(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event));
        Uerrors(k,event) = norm(relative_error,2)/sqrt(numel(relative_error));
        relative_error = (V(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event)-Vk(1:end-1,1:end-1,event))/abs(eps+V(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event));
        Verrors(k,event) = norm(relative_error,2)/sqrt(numel(relative_error));
        relative_error = (T(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event)-Tk(1:end-1,1:end-1,event))/abs(eps+T(1:2^(kf-k):end-1,1:2^(kf-k):end-1,event));
        Terrors(k,event) = norm(relative_error,2)/sqrt(numel(relative_error));
    end
end
%% plotting results
% load manually the .mat file with the data if variables not avaliable
close all;
kf = powers(end);
U = getfield(UU,strcat('k',num2str(kf)));
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
end 