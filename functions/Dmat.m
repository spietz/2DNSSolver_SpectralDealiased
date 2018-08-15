function [D,M] = Dmat(y,alpha,beta,N) 

%% generalised Vandemonde matrix
    
    % (transpose) Vandemonde matrix relating nodal to modal coefficietns
    V = zeros(N+1,N+1); % init matrix (i is modes, j is gl-nodes)
    %V = zeros(N-1,N-1); % init matrix (i is modes, j is gl-nodes)
    for n = 0:N % loop modes
        Pn = NormJacobiPol(y,alpha,beta,n); % normalized
        V(n+1,:) = Pn(:);
    end
    V = V';
    Vinv = inv(V);
    
    % (transpose) Vandemonde der. matrix relating nodal to modal coefficietns
    Vy = zeros(N+1,N+1); % init matrix (i is modes, j is gl-nodes)
    %Vy = zeros(N-1,N-1); % init matrix (i is modes, j is gl-nodes)
    for n = 0:N % loop modes
        dPndy = GradNormJacobiPol(y,alpha,beta,n,1); % normalized
        Vy(n+1,:) = dPndy(:);
    end
    Vy = Vy';
    
%% Differentiation matrix
    D = Vy*Vinv;
    
%% Mass matrix
    M = inv(V*V');