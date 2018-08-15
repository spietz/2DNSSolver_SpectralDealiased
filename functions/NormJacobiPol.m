function [P] = NormJacobiPol(x,alpha,beta,N);

% http://mathworld.wolfram.com/OrthogonalPolynomials.html
%Calculation of the constant
if ( alpha == -0.5 & beta == -0.5 )
    if N == 0
        norm_cnt = pi;
    else
        norm_cnt = 0.5*pi;
    end
elseif (N<95) 
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


P=JacobiPol(x,alpha,beta,N)/sqrt(norm_cnt);

