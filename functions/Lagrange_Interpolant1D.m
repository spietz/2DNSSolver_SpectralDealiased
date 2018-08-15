function INF=Lagrange_Interpolant1D(x,f,xgrid)

INF=0;
Nx=length(xgrid)-1; % is really Nx-2
for i=1:Nx+1
  carry = f(i)*Lagrange_poly(i-1,x,xgrid) ;
  INF=INF+carry;  
end
 
