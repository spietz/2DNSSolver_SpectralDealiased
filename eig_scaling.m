addpath('functions/')

% eigenvalues analysis

clc, clear all, close all

% Setup
xi              = 0.5;
eta             = 0.5;
Ulid            = 1;
beta            = sqrt(5.0);
Revec           = [1, 100, 400, 1000];
%Revec           = [2, 10, 25, 50, 100, 200, 400, 800, 1000, 2000, 5000];
Nvec            = [8, 16, 24, 32]; 

Nmax            = max(Nvec);       %max Nx
nBndmax         = 4*Nmax;          %max boundary points in grid
Nlambda_max     = 3*((Nmax-1)^2);  %max number of eigenvalues

LAMBDA  = zeros(length(Revec),length(Nvec),Nlambda_max);
LDT     = zeros(length(Revec),length(Nvec),Nlambda_max);
NLAMBDA = zeros(length(Revec),length(Nvec));

for r = 1:length(Revec)

    Re  = Revec(r)
    
    for n = 1:length(Nvec)
        
        Nx     = Nvec(n)
        Ny     = Nvec(n);
        
%         %Choosing max values for velocities
%         load(sprintf('../data/Re%i_Nx%i_Ny%i.mat',Re,Nx,Ny),...
%                'U','V')
%         close all
%         Umax           = max(max(U));
%         Vmax           = max(max(V));
        Umax            = Ulid; 
        Vmax            =    0;           
        
        [Rsmall,Rbig,idxBnd,X,Y]  = SpatialOperator(Nx, Ny, xi, eta,Umax,Vmax,Re);
              
        lambda          = eig(full(Rsmall));
              
        % storage eigenvalues
       nlambda          = length(lambda); %must be 3*(Nx-1)^2
       NLAMBDA(r,n)     = nlambda;          
       LAMBDA(r,n,1:length(lambda)) =lambda;
                            
    end
end

RL        = real(LAMBDA);
IL        = imag(LAMBDA);

maxRL     = max(abs(RL),[],3);
maxIL     = max(abs(IL),[],3);
% analisys of max eigen

%% Figure 1 .- maxRL scaling with Re 

fig1 = figure(1); clf
subplot(1,1,1)
clear str_list
for n=1:length(Nvec)
 loglog(Revec,maxRL(:,n),'.','Markersize',10)
 if(n==1), hold on, end
 str_list{n} = sprintf('$N_x =%i$',Nvec(n));
end
%auxiliar lines 
 loglog(Revec,500*Revec.^-1,'--k')    
 str_list{n+1}='$\mathcal{O}(R_e^{-1})$';
% loglog([Nvec 256],[Nvec 256].^-4,'-.k')    
% str_list{r+2}='$\mathcal{O}(N_x^{-4})$';
% xlim([0 300])

% use latex legends and labels, set fontsize/name
 legend(str_list,...
         'FontUnits','points',...
         'interpreter','latex',...
         'FontSize',14,...
         'FontName','Computer Modern Roman',...
         'Location','NorthWest');
  % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',14)
        
        xlabel('$R_e$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$max(Re(\lambda))$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
%print(fig1,'-depsc',sprintf('../Figs4Henrik/fig_scal_real_lambda_Rbig_Re'))

% Figure 2 .- maxIL scaling with Re 

fig2 = figure(2); clf
subplot(1,1,1)
clear str_list
for n=1:length(Nvec)
 loglog(Revec(1:end),maxIL(1:end,n),'.','Markersize',10)
 if(n==1), hold on, end
 str_list{n} = sprintf('$N_x =%i$',Nvec(n));
end
%auxiliar lines 
 loglog(Revec,50*Revec,'--k')    
 str_list{n+1}='$\mathcal{O}(Re)$';
 axis([1 2000 0 10^5])
% loglog([Nvec 256],[Nvec 256].^-4,'-.k')    
% str_list{r+2}='$\mathcal{O}(N_x^{-4})$';
% xlim([0 300])

% use latex legends and labels, set fontsize/name
 legend(str_list,...
         'FontUnits','points',...
         'interpreter','latex',...
         'FontSize',14,...
         'FontName','Computer Modern Roman',...
         'Location','SouthWest');
  % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$R_e$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$max(Im(\lambda))$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')

%print(fig2,'-depsc',sprintf('../Figs4Henrik/fig_scal_Imag_lambda_Rbig_Re'))
     
        
%% Figure 1 .- maxRL scaling with Nx

fig3 = figure(3); clf
subplot(1,1,1)
clear str_list
for r=1:length(Revec)
 loglog(Nvec,maxRL(r,:)','.','Markersize',15)
 if(r==1), hold on, end
 str_list{r} = sprintf('$Re =%i$',Revec(r));
end
%auxiliar lines 
 loglog(Nvec,Nvec.^4,'--k')    
 str_list{r+1}='$\mathcal{O}(N_x^{4})$';
% loglog([Nvec 256],[Nvec 256].^-4,'-.k')    
% str_list{r+2}='$\mathcal{O}(N_x^{-4})$';
% xlim([0 300])

% use latex legends and labels, set fontsize/name
 legend(str_list,...
         'FontUnits','points',...
         'interpreter','latex',...
         'FontSize',18,...
         'FontName','Computer Modern Roman',...
         'Location','NorthWest');
  % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$N_x$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$max(Re(\lambda))$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')        
 
 %       print(fig3,'-depsc',sprintf('../Figs4Henrik/fig_scal_real_lambda_Rbig_Nx'))
       
        
        
%% Figure 4 .- maxIL scaling with Nx

fig4 = figure(4); clf
subplot(1,1,1)
clear str_list
for r=1:length(Revec)
 if (r~=3)   
  loglog(Nvec,maxIL(r,:)','.','Markersize',15)
 else 
  loglog(Nvec,maxIL(r,:)','o','Markersize',12) 
 end 
 if(r==1), hold on, end
 str_list{r} = sprintf('$Re =%i$',Revec(r));
end
%auxiliar lines 
% loglog(Nvec,Nvec,'--k')    
% str_list{r+1}='$\mathcal{O}(N_x)$';
loglog(Nvec,5*Nvec.^2,'-.k')    
 str_list{r+1}='$\mathcal{O}(N_x^{2})$';
% xlim([0 300])

% use latex legends and labels, set fontsize/name
 legend(str_list,...
         'FontUnits','points',...
         'interpreter','latex',...
         'FontSize',18,...
         'FontName','Computer Modern Roman',...
         'Location','NorthEast');
  % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$N_x$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$max(Im(\lambda))$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')        
  
  %      print(fig4,'-depsc',sprintf('../Figs4Henrik/fig_scal_imag_lambda_Rbig_Nx'))
     
        