function [Te,U,V]=thermal_convection_projection()
%MIT18086_NAVIERSTOKES
%    Solves the incompressible Navier-Stokes equations in a
%    rectangular domain with prescribed velocities along the
%    boundary. The solution method is finite differencing on
%    a staggered grid with implicit diffusion and a Chorin
%    projection method for the pressure.
%    Visualization is done by a colormap-isoline plot for
%    pressure and normalized quiver and streamline plot for
%    the velocity field.
%    The standard setup solves a lid driven cavity problem.

% 07/2007 by Benjamin Seibold
%            http://www-math.mit.edu/~seibold/
% Feel free to modify for teaching and learning.
%-----------------------------------------------------------------------
dt = 3e-2;    % time step
nx = 65;      % number of x-gridpoints
ny = 17;      % number of y-gridpoints
nsteps = 1e6;  % number of steps with graphic output
%----------------------------- Thermal coefficients -------------------
%------------------------------Glycerine
% Lx = 0.38;                  % Cavity width [m]
% Ly = 0.04;                  % Cavity height [m]
% L = Ly;                     % Reference lenght [m]
% lx = Lx/L;                  % Normalized cavity width []
% ly = Ly/L;                  % Normalized cavity height []
% G = 9.8;                    % Gravity acceleration [m/s^2]
% g = 1*G;                    % Gravity acceleration [m/s^2]
% g = L*G;                    % Gravity acceleration [m/s^2]
% TC = 291.2;                 % Cold temperature [K]
% TH = 294.78;                % Hot temperature [K]
% T0 = 293;                   % Reference Temperature [K]
% Tc = (TC-T0)/(TH-TC);       % Cold temperature [K]
% Th = (TH-T0)/(TH-TC);       % Hot temperature [K]
% rho = 1264.02;              % Glycerine density [kg/m^3]
% mu = 1.49;                  % Dynamic viscosity []
% nu = mu/rho;                % Kinematic viscosity []
% B = 5e-4;                   % coeficient of thermal expansion [1/K]
% Pr = 1.25e4;                % Prandtl number
% u0 = 1;                     % Reference velocity [m/s]
% Re = rho*u0*L/mu;           % Reynolds number
% b = (Th-Tc)*B;              % normalized coeficient of thermal expansion []
% Gr = g*B*(Th-Tc)*(L^3/nu^2);% Grashof number
% Ra = Pr*Gr;                 % Rayleigh number
% tS = 1e4;                   % Real Experiment time [s]
% tf = tS/L;                  % Simulation time []
%------------------------------Air
Ly = 0.6557;                  % Cavity height [m]
Lx = Ly*4;                  % Cavity width [m]
L = Ly;                     % Reference lenght [m]
lx = Lx/L;                  % Normalized cavity width []
ly = Ly/L;                  % Normalized cavity height []
G = 9.8126;                    % Gravity acceleration [m/s^2]
%g = 1*G;                    % Gravity acceleration [m/s^2]
g = L*G;                   % Gravity acceleration [m/s^2]
TC = 292.5;                 % Cold temperature [K]
TH = 293.5;                % Hot temperature [K]
T0 = 293;                   % Reference Temperature [K]
Tc = (TC-T0)/(TH-TC);       % Cold temperature [K]
Th = (TH-T0)/(TH-TC);       % Hot temperature [K]
rho = 1.205;              % Glycerine density [kg/m^3]
mu = 1.81e-4;                  % Dynamic viscosity []
nu = mu/rho;                % Kinematic viscosity []
B = 3.4e-3;                   % coeficient of thermal expansion [1/K]
Pr = 0.72;                % Prandtl number
u0 = 1;                     % Reference velocity [m/s]
Re = rho*u0*L/mu;           % Reynolds number
b = (Th-Tc)*B;              % normalized coeficient of thermal expansion []
Gr = G*B*(Th-Tc)*(L^3/nu^2);% Grashof number
Ra = Pr*Gr;                 % Rayleigh number
tS = 1e4;                   % Real Experiment time [s] tS = 1e4; 
tf = tS/L;                  % Simulation time []

%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;
x = linspace(0,lx,nx+1); hx = lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
[X,Y] = meshgrid(y,x);
%-----------------------------------------------------------------------
% initial conditions
U = zeros(nx-1,ny); V = zeros(nx,ny-1); T = zeros(nx,ny-1)+Tc;
% boundary conditions
uN = x*0+0;    vN = avg(x)*0;   TN = avg(x)*0+Tc;
uS = x*0;      vS = avg(x)*0;   TS = avg(x)*0+Th;
uW = avg(y)*0; vW = y*0;        TW = y*0+0;
uE = avg(y)*0; vE = y*0;        TE = y*0+0;
%-----------------------------------------------------------------------
Ubc = dt/Re*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hx^2+...
      [uW;zeros(nx-3,ny);uE]/hy^2);
Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hx^2+...
      [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hy^2);
Tbc = ([TS' zeros(nx,ny-3) TN']/hx^2+...
      [2*TW(2:end-1);zeros(nx-2,ny-1);2*TE(2:end-1)]/hy^2);

fprintf('initialization')
Lp = kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx));
Lp(1,1) = 3/2*Lp(1,1);
perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';

Lu = speye((nx-1)*ny)+dt/Re*(kron(speye(ny),K1(nx-1,hx,2))+...
     kron(K1(ny,hy,3),speye(nx-1)));
peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
Lv = speye(nx*(ny-1))+dt/Re*(kron(speye(ny-1),K1(nx,hx,3))+...
     kron(K1(ny-1,hy,2),speye(nx)));
perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
LT = speye(nx*ny)+(1/(Re*Pr))*dt*(kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx)));
perT = symamd(LT); RT = chol(LT(perT,perT)); RTt = RT';
Lq = kron(speye(ny-1),K1(nx-1,hx,2))+kron(K1(ny-1,hy,2),speye(nx-1));
perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
figure(1);
fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')
for k = 1:nt
   % treat nonlinear terms
   gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
   Ue = [uW;U;uE]; Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
   Ve = [vS' V vN']; Ve = [2*vW-Ve(1,:);Ve;2*vE-Ve(end,:)];
   Te = [TS' T TN']; Te = [Te(1,:);Te;Te(end,:)];
   Ua = avg(Ue')'; Ud = diff(Ue')'/2;
   Va = avg(Ve);   Vd = diff(Ve)/2;
   UVx = diff(Ua.*Va-gamma*abs(Ua).*Vd)/hx;
   UVy = diff((Ua.*Va-gamma*Ud.*abs(Va))')'/hy;
   Ua = avg(Ue(:,2:end-1));   Ud = diff(Ue(:,2:end-1))/2;
   Va = avg(Ve(2:end-1,:)')'; Vd = diff(Ve(2:end-1,:)')'/2;
   U2x = diff(Ua.^2-gamma*abs(Ua).*Ud)/hx;
   V2y = diff((Va.^2-gamma*abs(Va).*Vd)')'/hy;
   U = U-dt*(UVy(2:end-1,:)+U2x);
   V = V-dt*(UVx(:,2:end-1)+V2y);
   % implicit viscosity
   rhs = reshape(U+Ubc,[],1);
   u(peru) = Ru\(Rut\rhs(peru));
   U = reshape(u,nx-1,ny);
   % Buoyancy force included
   rhs = reshape(V+Vbc+dt*(1-b*T)*g,[],1);
   %return

   v(perv) = Rv\(Rvt\rhs(perv));
   V = reshape(v,nx,ny-1);
   
   % pressure correction
   rhs = reshape(diff([uW;U;uE])/hx+diff([vS' V vN']')'/hy,[],1);
   p(perp) = -Rp\(Rpt\rhs(perp));
   P = reshape(p,nx,ny);
   U = U-diff(P)/hx;
   V = V-diff(P')'/hy;
   
   % temperature update
   %rhs = reshape(Ta - dt*Ua*diff([TW;T;TE])/hx+dt*Va*diff([TS' T TN']')'/hy,[],1);%Issue with Tx
   Tay = diff(Te(2:end-1,:)')';
   Ta = avg(Te(2:end-1,:)')';
   Tax = avg(avg(diff(Te))')';
   rhs = reshape(Ta - dt*Ua.*Tax/hx - dt*Va.*Tay/hy,[],1);%at cell centers
   t(perT) = RT\(RTt\rhs(perT));
   T = reshape(t,nx,ny);
   T = avg(T')';%Tshifted to match with v
   
   % visualization
   if floor(25*k/nt)>floor(25*(k-1)/nt), fprintf('.'), end
   if k==1|floor(nsteps*k/nt)>floor(nsteps*(k-1)/nt)
      % stream function
      rhs = reshape(diff(U')'/hy-diff(V)/hx,[],1);
      q(perq) = Rq\(Rqt\rhs(perq));
      Q = zeros(nx+1,ny+1);
      Q(2:end-1,2:end-1) = reshape(q,nx-1,ny-1);
      clf, contourf(avg(x),avg(y),Ta',20,'w-'), hold on
      contour(x,y,Q',20,'k-');
      Ue = [uS' avg([uW;U;uE]')' uN'];
      Ve = [vW;avg([vS' V vN']);vE];
      Len = sqrt(Ue.^2+Ve.^2+eps);
      quiver(x,y,(Ue./Len)',(Ve./Len)',.4,'k-');
      %quiver(x,y,Ue',Ve',.6,'k-');
      hold off, axis equal, axis([0 lx 0 ly])
      p = sort(p);% caxis(p([8 end-7]))
      title(sprintf('Re = %0.1g   t = %0.2g',Re,k*dt))
      colorbar;
      drawnow
   end
end
fprintf('\n')

%=======================================================================

function B = avg(A,k)
if nargin<2, k = 1; end
if size(A,1)==1, A = A'; end
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end
if size(A,2)==1, B = B'; end

function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;
