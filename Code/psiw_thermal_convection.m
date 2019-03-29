% psiw.m: Solves the unsteady driven cavity problem using the
%         streamfunction-vorticity algorithm described in
%         Pozrikidis (1997), pages 608-9.
%
% Author: John Stockie
%         Department of Mathematics
%         Simon Fraser University
%         stockie@math.sfu.ca
%
% Date:   February 11, 2019

% Physical constants:
Vd = 0;                 % lid velocity (m/s)
nsteps = 2e3;               % number of steps with graphic output
%------------------------------Glycerine
Lx = 0.38;                  % Cavity width [m]
Ly = 0.04;                  % Cavity height [m]
L = Ly;                     % Reference lenght [m]
lx = Lx/L;                  % Normalized cavity width []
ly = Ly/L;                  % Normalized cavity height []
G = 9.8;                    % Gravity acceleration [m/s^2]
g = 1*G;                    % Gravity acceleration [m/s^2]
%g = L*G;                   % Gravity acceleration [m/s^2]
TC = 291.2;                 % Cold temperature [K]
TH = 294.78;                % Hot temperature [K]
T0 = 293;                   % Reference Temperature [K]
Tc = (TC-T0)/(TH-TC);       % Cold temperature [K]
Th = (TH-T0)/(TH-TC);       % Hot temperature [K]
rho = 1264.02;              % Glycerine density [kg/m^3]
mu = 1.49;                  % Dynamic viscosity []
nu = mu/rho;                % Kinematic viscosity []
B = 5e-4;                   % coeficient of thermal expansion [1/K]
Pr = 1.25e4;                % Prandtl number
u0 = 1;                     % Reference velocity [m/s]
Re = rho*u0*L/mu;           % Reynolds number
b = (Th-Tc)*B;              % normalized coeficient of thermal expansion []
Gr = g*B*(Th-Tc)*(L^3/nu^2);% Grashof number
Ra = Pr*Gr;                 % Rayleigh number
tS = 1e4;                   % Real Experiment time [s] tS = 1e4; 
tend = tS/L;                  % Simulation time []
% Numerical parameters:
nx = 49;
ny = 5;

% Derived parameters:
dx = lx/nx;  
dy = ly/ny;

% Choose a time step based on the minimum of the two time step
% restrictions that come from the von Neumann analysis of the
% 2D linear advection-diffusion equation.
% dt1 = min(dx,dy) / Vd;              % advection restriction
% dt2 = 0.5 / nu / (1/dx^2 + 1/dy^2); % diffusion restriction
% dt3 = nu / Vd^2;                    % mixed (*not used*)
% safetyfac = 0.8;                    % "safety factor" (should be < 1)
% nt = floor(tend / (min(dt1,dt2) * safetyfac));
% dt = tend / nt;
dt = 3.2e-02;
nt = floor(tend / dt);

% Display parameters before starting.
% fprintf( 1, '\nReynolds number: %f\n', Re );
% fprintf( 1, 'Time step restrictions:\n    dt = %e (advective)\n', dt1 );
% fprintf( 1, '    dt = %e (diffusive)\n    dt = %e (mixed - not used)\n', dt2, dt3 );
% if dt1 < dt2, tstr = 'advection-limited';
% else          tstr = 'diffusion-limited'; end;
% fprintf( 1, 'Actual time step: %e  (%s)\n', dt, tstr );
%fprintf( 1, '\nPress <enter> to continue ...\n' );
%pause

% Set up an equally-spaced rectangular grid.
[xx,yy] = meshgrid(0:dx:lx, 0:dy:ly);

% STEP 1. Set up the initial conditions for psi, u and v, which are
% zero except for the lid velocity.  The unknowns are all
% node-centered and arrays are of size (nx+1)-by-(ny+1). 
psi = 0*xx; 
w   = psi;
u   = psi;
v   = psi;
T   = psi;
T(end,:) = Tc;
T(1,:) = Th;
u(end,:) = Vd; % lid BC

% Set up the matrix for the Poisson equation for psi.  The ordering
% of unknowns is row-wise: 
% (1,1), (1,2), (1,3), ... (1,Nx+1), (2,1), (2,2), ...
mx = 1/dx^2;
my = 1/dy^2;
e    = repmat( ones(nx-1,1), ny-1, 1);
esub = repmat( [ones(nx-2,1); 0], ny-1, 1);
esup = repmat( [0; ones(nx-2,1)], ny-1, 1);
nn = (nx-1)*(ny-1);
A = spdiags( [my*e, mx*esub, -2*(mx+my)*e, mx*esup, my*e], ...
	     [-(nx-1), -1, 0, 1, nx-1], nn, nn );

% These integer index vectors refer to interior grid points and are  
% useful when indexing arrays below.
ii = 2:nx;
jj = 2:ny;
fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')
for i = 1 : nt,
  % STEP 2. Compute vorticity by differencing the velocities,
  % first in the interior:
  w(jj,ii) = (v(jj,ii+1) - v(jj,ii-1)) / (2*dx) ...
      - (u(jj+1,ii) - u(jj-1,ii)) / (2*dy);
  % ... and then on the boundaries:
  w(jj,1)   = (4*v(jj,2) - v(jj,3)) / (2*dx);
  w(jj,end) = (v(jj,end-2) - 4*v(jj,end-1)) / (2*dx);
  w(1,ii)   = (-4*u(2,ii) + u(3,ii)) / (2*dy);
  w(end,ii) = (-u(end-2,ii) + 4*u(end-1,ii) - 3*Vd) / (2*dy);
  
  % STEP 3. Solve the vorticity equation using an explicit FTCS
  % scheme. 
  w(jj,ii) = w(jj,ii) - dt * ...
      ( u(jj,ii) .* (w(jj,ii+1) - w(jj,ii-1)) / (2*dx) ...
	+ v(jj,ii) .* (w(jj+1,ii) - w(jj-1,ii)) / (2*dy) ) ...
	+ (1/Re)*dt*(w(jj,ii+1) - 2*w(jj,ii) + w(jj,ii-1)) / dx^2 ...
	+ (1/Re)*dt*(w(jj+1,ii) - 2*w(jj,ii) + w(jj-1,ii)) / dy^2 ...
    + dt*b*(T(jj+1,ii) - T(jj-1,ii)) / (2*dy);
  
  % STEP 4. Solve the Poisson equation for streamfunction, by first 
  % setting up the right hand side for Delta psi = -w.  
  rhs = -w(jj,ii);  % this is a (nx-1) by (ny-1) matrix

  % Fix up the RHS for the boundary conditions.
  rhs(1,:)   = rhs(1,:)   - psi(1,ii)   / dy^2;
  rhs(end,:) = rhs(end,:) - psi(end,ii) / dy^2;
  rhs(:,1)   = rhs(:,1)   - psi(jj,1)   / dx^2;
  rhs(:,end) = rhs(:,end) - psi(jj,end) / dx^2;

  % Finally, convert the RHS to a vector and solve the linear system.
  rhs = reshape( rhs', nn, 1 ); % this is a (nx-1)*(ny-1)-vector
  psivec = A \ rhs;
  psi(jj,ii) = reshape(psivec, nx-1, ny-1)';

  % STEP 5. Compute the new velocity components by differencing
  % the streamfunction.  
  if floor(25*i/nt)>floor(25*(i-1)/nt), fprintf('.'), end
  u(jj,ii) =  (psi(jj+1,ii) - psi(jj-1,ii)) / (2*dy);
  v(jj,ii) = -(psi(jj,ii+1) - psi(jj,ii-1)) / (2*dx);
  
  % STEP 6.0 impose the boundary conditions boundaries for T:
  T(jj,1)   = T(jj,2);
  T(jj,end) = T(jj,end-1);
  T(1,ii)   = Th;
  T(end,ii) = Tc;
  
  % STEP 6.1 Compute the solution for the termal advection diffusion
  T(jj,ii) = T(jj,ii) - dt * ...
    ( u(jj,ii) .* (T(jj,ii+1) - T(jj,ii-1)) / (2*dx) ...
  + v(jj,ii) .* (T(jj+1,ii) - T(jj-1,ii)) / (2*dy) ) ...
  + (1/(Pr*Re))*dt*(T(jj,ii+1) - 2*T(jj,ii) + T(jj,ii-1)) / dx^2 ...
  + (1/(Pr*Re))*dt*(T(jj+1,ii) - 2*T(jj,ii) + T(jj-1,ii)) / dy^2;

%   if i==10
%   return
%   end

  % STEP 7. Present graphics
  if i==1|floor(nsteps*i/nt)>floor(nsteps*(i-1)/nt)
      figure(1);
      try 
          cla(plt1); cla(plt2);
      end
      plt1=subplot(1,2,1);streamslice(u,v);
      plt2=subplot(1,2,2);contourf(0:dx:lx,0:dy:ly,T);
      colorbar;
  end
end
