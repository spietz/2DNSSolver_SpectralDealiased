function [X, Y, GradX_vel, GradY_vel, Lap_vel, GradX_pres, GradY_pres, Mx_vel, My_vel] = operators(Nx, Ny, xi, eta)

%% Setup
if nargin < 5
    alpha  = 0;
    beta   = 0;
end

%% get nodes
xbar = JacobiGL(alpha,beta,Nx); % GL nodes
ybar = JacobiGL(alpha,beta,Ny);
% (u,v)-grid holds all points
% p-grid excludes end points

%% 2d grid
[X,Y] = meshgrid(xbar*xi,ybar*eta);

%% 1D differential operators
[Dx_vel, Mx_vel]  = Dmat(xbar,alpha,beta,Nx);
[Dy_vel, My_vel]  = Dmat(ybar,alpha,beta,Ny);
Dx_pres = zeros(Nx+1); % zero padding for consistency
Dx_pres(2:Nx,2:Nx) = Dmat(xbar(2:Nx),alpha,beta,Nx-2);
Dy_pres = zeros(Ny+1); % zero padding for consistency
Dy_pres(2:Ny,2:Ny) = Dmat(ybar(2:Ny),alpha,beta,Ny-2);

% scale
Dx_vel = Dx_vel/xi;
Dy_vel = Dy_vel/eta;
Dx_pres = Dx_pres/xi;
Dy_pres = Dy_pres/eta;
Mx_vel = Mx_vel*xi;
My_vel = My_vel*eta;

%% 2D differential operators
GradX_vel  = kron(Dx_vel,speye(Ny+1)); % x grad
GradY_vel  = kron(speye(Nx+1),Dy_vel); % y grad
Lap_vel    = GradX_vel*GradX_vel + GradY_vel*GradY_vel; % laplacian
Iy = sparse(Ny+1,Ny+1); Iy(2:Ny,2:Ny) = speye(Ny-1);
GradX_pres = kron(Dx_pres,Iy); % x grad
Ix = sparse(Nx+1,Nx+1); Ix(2:Nx,2:Nx) = speye(Nx-1);
GradY_pres = kron(Ix,Dy_pres); % y grad

GradX_vel = sparse(GradX_vel);
GradY_vel = sparse(GradY_vel);
Lap_vel = sparse(Lap_vel);
GradX_pres = sparse(GradX_pres);
GradY_pres = sparse(GradY_pres);
