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
Vd = 1.5;      % lid velocity (m/s)
nu = 0.1;      % kinematic viscosity (m^2/s)
Lx = 1.0;      % width of box (m)
Ly = 1.0;      % height of box (m)
tend = 0.2;    % length of time interval (s)
Re = Vd*Lx/nu; % Reynolds number (based on lid length)

% Numerical parameters:
nx = 40;
ny = 40;

% Derived parameters:
dx = Lx/nx;  
dy = Ly/ny;

% Choose a time step based on the minimum of the two time step
% restrictions that come from the von Neumann analysis of the
% 2D linear advection-diffusion equation.
dt1 = min(dx,dy) / Vd;              % advection restriction
dt2 = 0.5 / nu / (1/dx^2 + 1/dy^2); % diffusion restriction
dt3 = nu / Vd^2;                    % mixed (*not used*)
safetyfac = 0.8;                    % "safety factor" (should be < 1)
nt = floor(tend / (min(dt1,dt2) * safetyfac));
dt = tend / nt;

% Display parameters before starting.
fprintf( 1, '\nReynolds number: %f\n', Re );
fprintf( 1, 'Time step restrictions:\n    dt = %e (advective)\n', dt1 );
fprintf( 1, '    dt = %e (diffusive)\n    dt = %e (mixed - not used)\n', dt2, dt3 );
if dt1 < dt2, tstr = 'advection-limited';
else          tstr = 'diffusion-limited'; end;
fprintf( 1, 'Actual time step: %e  (%s)\n', dt, tstr );
fprintf( 1, '\nPress <enter> to continue ...\n' );
pause

% Set up an equally-spaced rectangular grid.
[xx,yy] = meshgrid(0:dx:Lx, 0:dy:Ly);

% STEP 1. Set up the initial conditions for psi, u and v, which are
% zero except for the lid velocity.  The unknowns are all
% node-centered and arrays are of size (nx+1)-by-(ny+1). 
psi = 0*xx; 
w   = psi;
u   = psi;
v   = psi;
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
	+ nu*dt*(w(jj,ii+1) - 2*w(jj,ii) + w(jj,ii-1)) / dx^2 ...
	+ nu*dt*(w(jj+1,ii) - 2*w(jj,ii) + w(jj-1,ii)) / dy^2;
  
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
  u(jj,ii) =  (psi(jj+1,ii) - psi(jj-1,ii)) / (2*dy);
  v(jj,ii) = -(psi(jj,ii+1) - psi(jj,ii-1)) / (2*dx);
  
  fprintf( 1, 't = %f\n', i*dt );

end
