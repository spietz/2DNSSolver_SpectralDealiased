addpath('functions/')

% Eigenvalues analysis

clc, clear all, close all

% Setup
xi              = 0.5;
eta             = 0.5;
Ulid            = 1;
beta            = sqrt(5.0);
Re              = 1;
Nx              = 16;
Ny              = Nx;

% % Choosing max values for velocities from the steady state
% % take for 
% Nx = 48;Ny=48;
% load(sprintf('../data/Re%i_Nx%i_Ny%i.mat',Re,Nx,Ny),...
%                'U','V')
% close all
% Umax             = max(max(U));
% Vmax             = max(max(V));
% %go back
% Nx = 8; Ny = Nx;

Umax            = Ulid; 
Vmax            =    0;    

% Spatial operator and eigenvalues for these Nx and Re

[Rsmall,Rbig,idxBnd,X,Y]  = SpatialOperator(Nx, Ny, xi, eta,Umax,Vmax,Re);
              
lambda          = eig(full(Rsmall));
rl              = real(lambda);
ri              = imag(lambda);

% HEURISTIC For the time step (conservatively, CFL<1)
dxmin           = min(diff(X(1,1:end)));
%dxmin           = min(diff(X(1,:)));
dymin           = min(diff(Y(1:end,1)));
%dymin           = min(diff(Y(:,1)));
lambda1         = (Umax + sqrt(Umax^2+beta^2))/dxmin + 1/(Re*dxmin^2);
lambda2         = (Vmax + sqrt(Vmax^2+beta^2))/dymin + 1/(Re*dymin^2);
dt              = 0.9/(lambda1+lambda2);
 
 
%% Figure of ERK4 linear stability region  

%axis
% x range 
x0              =  -3;
x1              = 0.5;
Nxx             = 301;
% y range 
y0              =  -3;
y1              =   3;
Nyy             = 301;
% mesh
xv              = linspace(x0,x1,Nxx);
yv              = linspace(y0,y1,Nyy);
[x,y]           = meshgrid(xv,yv);
% Calculate z
z               = x+1i*y;
% 4th order Runge-Kutta growth factor
g4              = 1+z+0.5*z.^2+1/6*z.^3+1/24*z.^4;
% magnitude of g4
g4mag           = abs(g4);

figure(1)
contour(x,y,g4mag,[1 1],'k-');
set(gca,'FontSize',18)
axis([x0,x1,y0,y1]);
axis('square');
xlabel('Real $\lambda_i\Delta t$',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',18,...
        'FontName','Times');
ylabel('Imag $\lambda_i\Delta t$',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',18,...
        'FontName','Times');
grid on;
%print(gcf,'-depsc','RK4.eps');



 hold on

plot(rl*dt,ri*dt,'d') %umax in red
