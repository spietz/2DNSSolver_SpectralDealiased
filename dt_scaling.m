addpath('functions/')

% dt analysis

clc, clear all, close all

% Setup
xi              = 0.5;
eta             = 0.5;
Ulid            = 1;
beta            = sqrt(5.0);
Revec   = [1, 25, 100, 400, 1000];
Nvec    = [8, 16, 32, 48, 56]; 

DT      = zeros(length(Revec),length(Nvec));

for r = 1:length(Revec)

    Re  = Revec(r)
    
    for n = 1:length(Nvec)
        
        Nx     = Nvec(n)
        Ny     = Nvec(n);
        
%         %Choosing max values for velocities
%         load(sprintf('../data/Re%i_Nx%i_Ny%i.mat',Re,Nx,Ny),...
%                'U','V')
%         close all
%         Umax   = max(max(U));
%         Vmax   = max(max(V));
        Umax    = Ulid; 
        Vmax    =    0;           
        
        xbar      = JacobiGL(0,0,Nx);
        x         = xbar*xi;
        ybar      = JacobiGL(0,0,Ny);
        y         = ybar*eta;      
            
        % Choose time step (conservatively, CFL<1)
        dxmin           = min(diff(x));
        dymin           = min(diff(y));
        lambda1         = (Umax + sqrt(Umax^2+beta^2))/dxmin + 1/(Re*dxmin^2);
        lambda2         = (Vmax + sqrt(Vmax^2+beta^2))/dymin + 1/(Re*dymin^2);
        dt              = 0.9/(lambda1+lambda2);
            
        %Storage the steps
        DT(r,n)          = dt;                   
    end
end

%% Figure 1 .- dt scaling with Nx 

fig2 = figure(2); clf
subplot(1,1,1)
for r=1:length(Revec)
 if (r==4)
  loglog(Nvec,DT(r,:),'o','Markersize',10)
 else   
  loglog(Nvec,DT(r,:),'.','Markersize',10)
 end
 if(r==1), hold on, end
 str_list{r} = sprintf('Re =%i',Revec(r));
end
%auxiliar lines 
loglog([Nvec 256],[Nvec 256].^-2,'--k')    
str_list{r+1}='$\mathcal{O}(N_x^{-2})$';
loglog([Nvec 256],[Nvec 256].^-4,'-.k')    
str_list{r+2}='$\mathcal{O}(N_x^{-4})$';
axis([0 100 1.e-08 10])

% use latex legends and labels, set fontsize/name


 legend(str_list,...
         'FontUnits','points',...
         'interpreter','latex',...
         'FontSize',12,...
         'FontName','Computer Modern Roman',...
         'Location','NorthWest');
  % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',12)
        
        xlabel('$N_x$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$\Delta \tau$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
