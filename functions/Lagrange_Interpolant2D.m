function INF=Lagrange_Interpolant2D(x,y,F,xgrid,ygrid)

INF=zeros(length(y),length(x));
Npx=length(xgrid)-1; % is really Nx-2
Npy=length(ygrid)-1;
for i=1:Npx+1
 for j=1:Npy+1   
   carry = F(j,i)*Lagrange_poly(j-1,y,ygrid)*Lagrange_poly(i-1,x,xgrid)';
   INF=INF+carry;
 end
end

