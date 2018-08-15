function [P] = GradNormJacobiPol(x,alpha,beta,N,k)
% slides " Lecture 3, Polynomial Methods" p. 20
% derivatives of orthonormal JacobiPol!

if N >= k
    P = gamma(alpha+beta+N+1+k)/ (2^k*gamma(alpha+beta+N+1) ) ...
        * sqrt( norm_constant(alpha+k, beta+k, N-k) ...
        / norm_constant(alpha, beta, N) ) ...
        * NormJacobiPol(x,alpha+k,beta+k,N-k);
else
    P = zeros(size(x)); % derivative of constant is 0
end

end

function norm_cnt = norm_constant(alpha,beta,N)
%Calculation of the constant
if (N<95)
    norm_cnt=2^(alpha+beta+1) * (gamma(N+alpha+1)*gamma(N+beta+1)) / ...
        (factorial(N)*(2*N+alpha+beta+1)*gamma(N+alpha+beta+1));
else %we need to solve a limit, valid only for alpha beta integer
    prod1=1;prod2=1;
    for i=1:alpha
        prod1=prod1*(N+i);
        prod2=prod2*(N+beta+i);
    end
    norm_cnt=(2^(alpha+beta+1) / (2*N+alpha+beta+1)) * ...
        prod1/prod2;
end
end