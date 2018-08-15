addpath('functions/')

% driver

clc, clear all, close all

% Setup
dealias = true
restart = false
xi = 0.5;
eta = 0.5;
Ulid = 1;
beta = sqrt(5.0);
Revec = [1000]%, 100, 400, 1000];
Nvec = [8];
maxstep = 1000000;
steadytol = 1e-6;
statstride = 50;
step0 = 1;

for r = 1:length(Revec)
    Re=Revec(r)
    
    for nv = 1:length(Nvec)
        Nx = Nvec(nv)
        Ny = Nvec(nv);
        
        % Choose start guess
        P = zeros(Ny+1,Nx+1); dP = zeros(Ny+1,Nx+1);
        U = zeros(Ny+1,Nx+1); dU = zeros(Ny+1,Nx+1);
        V = zeros(Ny+1,Nx+1); dV = zeros(Ny+1,Nx+1);
        
        if restart
            % load selected variables
            load(sprintf('data/Re%i_Nx%i_Ny%i.mat',Re,Nx,Ny),...
                'U','V','P','step','steadyhist','statstride')
            Umax = max(abs(U(:)));
            Vmax = max(abs(V(:)));
            step0 = step + 1;
            steadyhist = [steadyhist, zeros(2,floor((maxstep-step0+1)/statstride),1)];
            
        else
            
            % smooth bc
            U(Ny+1,3:Nx-1) = Ulid; U(Ny+1,2) = Ulid/2; U(Ny+1,Nx) = Ulid/2;
            % U(Ny+1,1:Nx+1) = Ulid; [PSF, center] = psfGauss([Ny+1,Nx+1],2);
            % U = conv2(U,PSF,'same');
            
            Umax = Ulid;
            Vmax = 0;
            steadyhist = zeros(2,maxstep/statstride,1);
        end
        
        tstart = tic;
        
        % Generate grid and operators
        [X, Y, GradX_vel, GradY_vel, Lap_vel, GradX_pres, GradY_pres, Mx, My] ...
            = operators(Nx, Ny, xi, eta);
        N = (Nx+1)*(Ny+1);
        dxmin = min(diff(X(1,:)));
        dymin = min(diff(Y(:,1)));
        
        %% DEALIASING
        
        % rescale coarse axis to [-1,1]
        x1 = X(1,:)/xi; y1 = Y(:,1)/eta;
        
        % generate fine grid axis
        [x2, wx2] = JacobiGQ(0,0,2*Nx);
        [y2, wy2] = JacobiGQ(0,0,2*Ny);
        [X2,Y2]   = meshgrid(x2,y2);
        
        % coarse to fine pol. interpolation matrix (kopriva algorithm)
        % 1D
        Imatx = PolynomialInterpolationMatrix(x1,x2);
        Imaty = PolynomialInterpolationMatrix(y1,y2);
        % 2D
        ImatXY = kron(Imatx,Imaty);
        
        % allocate
        U2 = zeros(2*Ny+1, 2*Nx+1);
        V2 = zeros(2*Ny+1, 2*Nx+1);
        
        Uprod2 = zeros(2*Ny+1, 2*Nx+1);
        Vprod2 = zeros(2*Ny+1, 2*Nx+1);
        
        Uprodtx2 = zeros(2*Ny+1, 2*Nx+1);
        Vprodtx2 = zeros(2*Ny+1, 2*Nx+1);
        
        Uprodtxy2 = zeros(2*Ny+1, 2*Nx+1);
        Vprodtxy2 = zeros(2*Ny+1, 2*Nx+1);
        
        Uprodtx1 = zeros(Ny+1, Nx+1);
        Vprodtx1 = zeros(Ny+1, Nx+1);
        
        Uprod1 = zeros(Ny+1, Nx+1);
        Vprod1 = zeros(Ny+1, Nx+1);
        
        %% Pseudo time stepping
        w = [1/4,1/3,1/2,1]; % RK weights
        for step = step0:maxstep
            
            % Choose time step (conservatively)
            lambda1 = (Umax + sqrt(Umax^2+beta^2))/dxmin + 1/(Re*dxmin^2);
            lambda2 = (Vmax + sqrt(Vmax^2+beta^2))/dymin + 1/(Re*dymin^2);
            dt = 0.9/(lambda1+lambda2);
            
            Pold = P;
            Uold = U;
            Vold = V;
            
            for s = 1:4
                
                if dealias
                    %% Compute dealiased convection (Uprod1, Vprod1)
                    % interp vel and grad_vel to fine grid and compute product
                    U2(:)  = ImatXY*U(:);
                    V2(:)  = ImatXY*V(:);
                    Uprod2(:) = - U2(:).*( ImatXY*(GradX_vel*U(:)) ) ...
                        - V2(:).*( ImatXY*(GradY_vel*U(:)) );
                    Vprod2(:) = - U2(:).*( ImatXY*(GradX_vel*V(:)) ) ...
                        - V2(:).*( ImatXY*(GradY_vel*V(:)) );
                    
                    % discrete polynomial transform
                    % transforms first direction
                    for n = 0:2*Nx
                        nn = n + 1;
                        Pn = NormJacobiPol(x2,0,0,n); % column
                        Uprodtx2(:, nn) = Uprod2*(Pn.*wx2);
                        Vprodtx2(:, nn) = Vprod2*(Pn.*wx2);
                    end
                    % transform other direction (tranpose)
                    for m = 0:2*Ny
                        mm = m + 1;
                        Pm = NormJacobiPol(y2,0,0,m); % column
                        Uprodtxy2(mm, :) = (Uprodtx2'*(Pm.*wy2))';
                        Vprodtxy2(mm, :) = (Vprodtx2'*(Pm.*wy2))';
                    end
                    
                    % projection (truncated polynomial evaluated on coarse grid)
                    % transforms first direction
                    Uprodtx1(:,:) = 0; Vprodtx1(:,:) = 0;
                    for n = 0:Nx
                        nn = n + 1;
                        Pn = NormJacobiPol(x1,0,0,n); % column
                        Uprodtx1 = Uprodtx1 + Uprodtxy2(1:Ny+1, nn)*Pn';
                        Vprodtx1 = Vprodtx1 + Vprodtxy2(1:Ny+1, nn)*Pn';
                    end
                    % transform other direction (tranpose)
                    Uprod1(:,:) = 0; Vprod1(:,:) = 0;
                    for m = 0:Ny
                        mm = m + 1;
                        Pm = NormJacobiPol(y1,0,0,m); % column
                        Uprod1 = Uprod1 + (Uprodtx1(mm, 1:Nx+1)'*Pm')';
                        Vprod1 = Vprod1 + (Vprodtx1(mm, 1:Nx+1)'*Pm')';
                    end
                else
                    Uprod1 = - U(:).*(GradX_vel*U(:)) - V(:).*(GradY_vel*U(:));
                    Vprod1 = - U(:).*(GradX_vel*V(:)) - V(:).*(GradY_vel*V(:));
                end
                %% pseudo time stage derivatives
                dP(:) = -(beta^2)*GradX_vel*U(:) ...
                    -(beta^2)*GradY_vel*V(:);
                dU(:) = -GradX_pres*P(:)...
                    + (1/Re*Lap_vel)*U(:)...
                    + Uprod1(:);
                dV(:) = -GradY_pres*P(:)...
                    + (1/Re*Lap_vel)*V(:)...
                    + Vprod1(:);
                
                % advance pressure
                P = Pold + w(s)*dt*dP;
                
                % advance momentum
                U(2:Ny,2:Nx) = Uold(2:Ny,2:Nx) + w(s)*dt*dU(2:Ny,2:Nx);
                V(2:Ny,2:Nx) = Vold(2:Ny,2:Nx) + w(s)*dt*dV(2:Ny,2:Nx);
                
            end
            
            dU = (U - Uold);
            dV = (V - Vold);
            
            changeU = norm(dU(:),2)/norm(Uold(:),2);
            changeV = norm(dV(:),2)/norm(Vold(:),2);
            
            Umax = max(abs(U(:)));
            Vmax = max(abs(V(:)));
            
            if max(changeU,changeV) < steadytol
                break
            end
            
            if mod(step, statstride) == 0
                step
                changeU
                changeV
                dt;
                steadyhist(1,step/statstride) = changeU;
                steadyhist(2,step/statstride) = changeV;
            end
            
        end
        
        telapsed = toc(tstart);
        
        %% Post processing
        
        % calc vorticity
        VORT = zeros(Ny+1,Nx+1);
        VORT(:) = GradX_vel*V(:) - GradY_vel*U(:);
        
        DVORTx = zeros(Ny+1,Nx+1);
        DVORTy = zeros(Ny+1,Nx+1);
        DVORTx(:) = GradX_vel*VORT(:);
        DVORTy(:) = GradY_vel*VORT(:);
        
        % Kinetic Energy
        E = 0.5*integrate2d(U.^2,Nx,Ny,Mx,My)
        
        % Enstrophy
        Z = 0.5*integrate2d(VORT.^2,Nx,Ny,Mx,My)
        
        % Palinenstrophy
        Pa = 0.5*integrate2d(DVORTx.^2+DVORTy.^2,Nx,Ny,Mx,My)
        
        % Save entire work space
        save(sprintf('data/Dealiased%i_Re%i_Nx%i_Ny%i.mat',dealias,Re,Nx,Ny))
        
        %% plots
        close all
        
        fig1 = figure(1); clf
        subplot(1,1,1)
        surf(X,Y,U)
        xlabel('x')
        ylabel('y')
        view([0 90])
        colorbar
        axis equal
        
        % % use latex legends and labels, set fontsize/name
        % legend({'$y=x$'},...
        %        'FontUnits','points',...
        %        'interpreter','latex',...
        %        'FontSize',18,...
        %        'FontName','Computer Modern Roman',...
        %        'Location','EO');
        % % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$x$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$y$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        print(fig1,'-depsc',sprintf('fig/Dealiased%i_U_Re%i_Nx%i_Ny%i',dealias,Re,Nx,Ny))
        
        fig2 = figure(2); clf
        subplot(1,1,1)
        surf(X,Y,V)
        xlabel('x')
        ylabel('y')
        view([0 90])
        colorbar
        axis equal
        
        % % use latex legends and labels, set fontsize/name
        % legend({'$y=x$'},...
        %        'FontUnits','points',...
        %        'interpreter','latex',...
        %        'FontSize',18,...
        %        'FontName','Computer Modern Roman',...
        %        'Location','EO');
        % % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$x$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$y$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        print(fig2,'-depsc',sprintf('fig/Dealiased%i_V_Re%i_Nx%i_Ny%i',dealias,Re,Nx,Ny))
        
        fig3 = figure(3); clf
        subplot(1,1,1)
        surf(X(2:Ny,2:Nx),Y(2:Ny,2:Nx),P(2:Ny,2:Nx))
        xlabel('x')
        ylabel('y')
        view([0 90])
        colorbar
        axis equal
        
        % % use latex legends and labels, set fontsize/name
        % legend({'$y=x$'},...
        %        'FontUnits','points',...
        %        'interpreter','latex',...
        %        'FontSize',18,...
        %        'FontName','Computer Modern Roman',...
        %        'Location','EO');
        % % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$x$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$y$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        print(fig3,'-depsc',sprintf('fig/Dealiased%i_P_Re%i_Nx%i_Ny%i',dealias,Re,Nx,Ny))
        
        fig4 = figure(4); clf
        subplot(1,1,1)
        surf(X,Y,VORT)
        xlabel('x')
        ylabel('y')
        view([0 90])
        colorbar
        axis equal
        
        % % use latex legends and labels, set fontsize/name
        % legend({'$y=x$'},...
        %        'FontUnits','points',...
        %        'interpreter','latex',...
        %        'FontSize',18,...
        %        'FontName','Computer Modern Roman',...
        %        'Location','EO');
        % % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$x$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$y$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        print(fig4,'-depsc',sprintf('fig/Dealiased%i_VORT_Re%i_Nx%i_Ny%i',dealias,Re,Nx,Ny))
        
        fig5 = figure(5); clf
        loglog((1:length(steadyhist))*statstride, steadyhist)
        
        % use latex legends and labels, set fontsize/name
        legend({'$u$', '$v$'},...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Computer Modern Roman',...
            'Location','EO');
        % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$n$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$R$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        print(fig5,'-depsc',sprintf('fig/Dealiased%i_Change_Re%i_Nx%i_Ny%i',dealias,Re,Nx,Ny))
        
        fig6 = figure(6); clf
        streamslice(X,Y,U,V), axis tight
        
        % use latex legends and labels, set fontsize/name
        legend({'$u$', '$v$'},...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Computer Modern Roman',...
            'Location','EO');
        % legend boxoff % if not outside plot?
        
        set(gca,'FontSize',18)
        
        xlabel('$x$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        ylabel('$y$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',18,...
            'FontName','Times')
        
        print(fig6,'-depsc',sprintf('fig/Dealiased%i_STREAMSLICE_Re%i_Nx%i_Ny%i',dealias,Re,Nx,Ny))
        
        %%
    end
end