addpath('functions/')

% driver

clc, clear all, close all

% setup
Re = 1000
N = 16;

% load simulation
load(sprintf('data/dealiased/Dealiased1_Re%i_Nx%i_Ny%i.mat',Re,N,N))
close all % close the figures that were saved

disp('******************************************************************');
disp(...
    sprintf('Loaded Re=%i, (Nx,Ny)=(%i,%i)', Re, Nx, Ny)...
    );
disp(...
    sprintf('final step/maxstep %i/%i', step, maxstep)...
    );
disp(...
    sprintf('with last saved residual (u,v) was (%3.2e, %3.2e)',steadyhist(1,floor(step/statstride)),steadyhist(2,floor(step/statstride)))...
    );

%% calc vorticity
VORT = zeros(Ny+1,Nx+1);
VORT(:) = GradX_vel*V(:) - GradY_vel*U(:);

%% compute streamfunction
xbar=JacobiGL(0,0,Nx);
ybar=JacobiGL(0,0,Ny);

% Interpolation centerlines grids
% we interpolate in
points=128;
xi=linspace(-1,1,points)'; % Working interval [-1,1]
yi=xi;
[Xi, Yi] = meshgrid(xi, yi);

% Ui = zeros(size(Xi)); Vi = Ui;
Ui = Lagrange_Interpolant2D(xi,yi,U,xbar,ybar);
Vi = Lagrange_Interpolant2D(xi,yi,V,xbar,ybar);

psi = cumsum(Ui*1/points);

%% plots

psilevels = [ -0.1, -0.08, -0.06, -0.04, -0.02, -0.01, -3e-3, -1e-3, -3e-4, -1e-4, -3e-5, -1e-5, -3e-6, -1e-6, -1e-7, -1e-8, -1e-9, -1e-10, 0.0  , 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01 , 0.03 , 0.05 , 0.07 , 0.09, 0.1  , 0.11 , 0.115, 0.1175];
vortlevels = [-40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.0, -0.5, -0.2, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]; 

fig1 = figure(1), clf

subplot(1,2,1)
hold on
contour(Xi*0.5,Yi*0.5,psi,psilevels,'k') % vorticity
hold off
axis off, axis equal, axis tight

subplot(1,2,2)
hold on
colorDepth = 5;
colormap(b2r(-colorDepth,colorDepth));
pcolor(X*0.5,Y*0.5,VORT*0.5) % vorticity
shading flat; % do not interpolate pixels
contour(X*0.5,Y*0.5,VORT,vortlevels,'k') % vorticity
hold off
axis off, axis equal, axis tight

print(fig1,'-depsc',sprintf('fig/Deliased1_stream_vort_Re%i_Nx%i_Ny%i',Re,Nx,Ny))