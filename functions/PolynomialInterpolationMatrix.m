function Imat = PolynomialInterpolationMatrix(X,Xq)
% Kopriva (and I):
% Algorithm 32:
% 1D PolynomialInterpolationMatrix:
% Matrix for Interpolation Between Two Sets of Points
% Uses Algorithms:
% Algorithm 30 (BarycentricWeights)
% Algorithm 139 (AlmostEqual)

X = X(:); % from-grid
Xq = Xq(:); % to-grid
N = length(X); % size of from-grid
M = length(Xq); % size of to-grid
Imat = zeros(M,N); % interpolation matrix
W = BarycentricWeights(X); % Weights for Lagrange Interpolation

for k = 1:M
    rowHasMatch = false;
    for j = 1:N
        Imat(k,j) = 0;
        if AlmostEqual(Xq(k),X(j))
            rowHasMatch = true;
            Imat(k,j) = 1;
        end
    end
    if rowHasMatch == false
        s = 0;
        for j = 1:N
            t = W(j)/(Xq(k)-X(j));
            Imat(k,j) = t;
            s = s + t;
        end
        for j = 1:N
            Imat(k,j) = Imat(k,j)/s;
        end
    end
end
Imat = sparse(Imat); % quick and dirty

function w = BarycentricWeights(x)
% Kopriva (and I):
% Algorithm 30:
% BarycentricWeights:
% Weights for Lagrange Interpolation

N = length(x);

for j = 1:N
    w(j) = 1;
end
for j = 2:N
    for k = 1:j-1
        w(k) = w(k)*(x(k)-x(j));
        w(j) = w(j)*(x(j)-x(k));
    end
end
for j = 1:N
    w(j) = 1/w(j);
end

function ae = AlmostEqual(a, b, tol)
% Kopriva (and I): 
% Algorithm 139: 
% AlmostEqual: Testing Equality of Two Floating Point Numbers
if nargin < 3
    tol = sqrt(eps);
end
if a == 0 || b == 0
    if abs(a-b) <= 2*tol
        ae = true;
    else
        ae = false;
    end
else
    if abs(a-b) <= abs(a)*tol && abs(a-b) <= abs(b)*tol
        ae = true;
    else
        ae = false;
    end
end
