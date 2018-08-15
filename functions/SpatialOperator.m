function [Rsmall,Rbig,idxBnd,X,Y] = SpatialOperator(Nx, Ny, xi, eta,Umax,Vmax,Re)

%% Setup
alpha = 0; beta = 0; % Legendre

%% get nodes
xbar = JacobiGL(alpha,beta,Nx); % GL nodes
ybar = JacobiGL(alpha,beta,Ny);
% (u,v)-grid holds all points
% p-grid excludes end points

%% 2d grid
[X,Y] = meshgrid(xbar*xi,ybar*eta);


%% Build index maps for grid points
M=Nx;N=Ny;
MN          = (M+1)*(N+1); % total # grid points
idx         = zeros((N+1),(M+1)); % initialize idx-mat
idx(:)      = 1:MN; % grid point #-mat
%  West and East boundaries are located on idxW, idxE
idxW        = idx(:,1);          idxW = idxW(:); 
idxE        = idx(:,end);        idxE = idxE(:); 
%  South and North boundaries on idxS, idxN
idxS        = idx(1,:) ;         idxS = idxS(:); 
idxN        = idx(end,:) ;       idxN = idxN(:);  
%  We stack them together in idxBnd ordered NSEW
idxBnd = [idxN;idxS;idxE;idxW]; %the corners are repeated here
idxBnd=unique(idxBnd);          %now idxBnd is sort and with no repetitions

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

%% the Big spatial operator R needs frozen aproaches Umax and Vmax 

diff_conv = - Umax*GradX_vel - Vmax*GradY_vel + 1/Re*Lap_vel;
beta = sqrt(5);
Rpp = zeros(MN); Ruv = zeros(MN); Rvu = zeros(MN);
Rpu         = -(beta^2)*GradX_vel;
Rpv         = -(beta^2)*GradY_vel;
Rup         = -GradX_pres;
Ruu         = diff_conv;
Rvp         = -GradY_pres;
Rvv         = diff_conv;
%constructing the matrix
R           = [Rpp Rpu Rpv ;
               Rup Ruu Ruv ;   % 3N x 3N matrix
               Rvp Rvu Rvv ]; 
            
idxchange=[idxBnd;MN+idxBnd;2*MN+idxBnd]; 

%Rsmall we erase those rows and columns regarding Boundary points
Rsmall = R;

Rsmall(idxchange,:)=[]; % we eliminate those rows
Rsmall(:,idxchange)=[]; % and columns  
           
Rsmall=sparse(Rsmall);

%Rbig we keep the size of R but replace the rows corrsp to indexBoundary to [0 0 ... 1 ...0 0]
Rbig = R;

Rbig(idxchange,:)            = zeros(length(idxchange), 3*MN);
Rbig(:,idxchange)            = zeros(3*MN,length(idxchange));

for i=1:length(idxchange)
 this_one = idxchange(i);   
 Rbig(this_one,this_one)    = 1 ;
end

Rbig                         = sparse(Rbig);

