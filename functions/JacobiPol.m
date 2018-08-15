function [P] = JacobiPol(x,alpha,beta,N);

% always return column
x = x(:);

% Jacobi polynomials for x in [-1,1] by recursion
P = zeros(length(x),3);

% starting values
P(:,1) = 1; % n = 0
P(:,2) = 0.5*(alpha-beta+(alpha+beta+2)*x); % n = 1

% calculate Pn+1
for n = 1:N-1
    
    a1 = 2*(n+alpha)*(n+beta)/((2*n+alpha+beta+1)*(2*n+alpha+beta)); % n-1,n
    a2 = (alpha^2-beta^2)/((2*n+alpha+beta+2)*(2*n+alpha+beta)); % n,n
    a3 = 2*(n+1)*(n+alpha+beta+1)/((2*n+alpha+beta+2)*(2*n+alpha+beta+1)); % n+1,n
    
    P(:,3) = ( ( a2 + x ) .* P(:,2) - a1 * P(:,1) ) / a3; % Pn+1
    
    % update
    P(:,1) = P(:,2); % n
    P(:,2) = P(:,3); % n - 1
end

if N < 2
    P = P(:,N+1);
else
    P = P(:,3);
end