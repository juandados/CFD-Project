% Natural convective heat transfer solver
% in a square domain
clc,
clear,
clf,

% Physical constants:
T0 = 0;                % Low wall temperature
Tw = 1;                % High  wall temperature 
Lx = 10.0;             % width of box (m)
Ly = 1.0;              % height of box (m)
tend = 5;              % length of time interval (s)

Pr = 0.7;              %Prandtl number
Ra = 2000;             %Rayleigh number

 
 

% Numerical parameters:
nx = 20;
ny = 20;

% Derived parameters:
dx = Lx/nx;  
dy = Ly/ny;

dt = 10^-4;
nt = tend/dt;

% Set up an equally-spaced rectangular grid.
[xx,yy] = meshgrid(0:dx:Lx, 0:dy:Ly);

% STEP 1. Set up the initial conditions for psi, T, u and v,
% All sit at node-centered and arrays are of size (nx+1)-by-(ny+1). 
psi = 0*xx; 
w   = psi;
u   = psi;
v   = psi;
T   = T0*ones(nx+1,ny+1);
T(1,:) = Tw;  % High temperature wall

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

for i = 1:nt
  
  t = i*dt;
  
  % STEP 2. Compute vorticity by differencing the velocities,
  % first in the interior:
  w(jj,ii) = (v(jj,ii+1) - v(jj,ii-1)) / (2*dx) ...
      - (u(jj+1,ii) - u(jj-1,ii)) / (2*dy);
  % ... and then on the boundaries:
  w(jj,1)   = (4*v(jj,2) - v(jj,3)) / (2*dx);
  w(jj,end) = (v(jj,end-2) - 4*v(jj,end-1)) / (2*dx);
  w(1,ii)   = (-4*u(2,ii) + u(3,ii)) / (2*dy);
  w(end,ii) = (-u(end-2,ii) + 4*u(end-1,ii)) / (2*dy);
  
  % STEP 3a. Solve the vorticity equation using an explicit FTCS
  % scheme. 
  w(jj,ii) = w(jj,ii) - dt * ...
      ( u(jj,ii) .* (w(jj,ii+1) - w(jj,ii-1)) / (2*dx) ...
	+ v(jj,ii) .* (w(jj+1,ii) - w(jj-1,ii)) / (2*dy) ) ...
	+ (Pr)*dt*(w(jj,ii+1) - 2*w(jj,ii) + w(jj,ii-1)) / dx^2 ...
	+ (Pr)*dt*(w(jj+1,ii) - 2*w(jj,ii) + w(jj-1,ii)) / dy^2 ...
    + (Pr*Ra)*(T(jj,ii+1) - T(jj,ii-1)) / dx;

  % STEP 3b. Solve the energy equation using an explicit FTCS
  T(jj,ii) = T(jj,ii) - dt * ...
      ( u(jj,ii) .* (T(jj,ii+1) - T(jj,ii-1)) / (2*dx) ...
	+ v(jj,ii) .* (T(jj+1,ii) - T(jj-1,ii)) / (2*dy) ) ...
	+ (1)*dt*(T(jj,ii+1) - 2*T(jj,ii) + T(jj,ii-1)) / dx^2 ...
	+ (1)*dt*(T(jj+1,ii) - 2*T(jj,ii) + T(jj-1,ii)) / dy^2; ...  
  
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
  
  figure(1)
  pcolor(xx,yy,T)
  xlabel('x'), ylabel('y')
  shading interp
  colorbar
  title(['Temperature, t=',num2str(t,'%f')])
  
  figure(2)
  quiver(xx,yy,u,v,'AutoScaleFactor',3,'MaxHeadSize',1)
  xlabel('x'), ylabel('y')
  xlim([0 Lx]);ylim([0 Ly]);   
  title(['Velocity, t=',num2str(t,'%f')])
  
  figure(3)
  contour(xx,yy,psi)
  xlabel('x'), ylabel('y')
  xlim([0 Lx]);ylim([0 Ly]);   
  title(['StreamLine, t=',num2str(t,'%f')])
  

end 

