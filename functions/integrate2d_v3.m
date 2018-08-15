function [fintxy] = integrate2d_v3(f, Nx, Ny, Mx, My)
% remember we always integrate |integrant|² using discrete matrix

%Mask construction
fplus=zeros(Ny+1,Nx+1); fminus=fplus;
for i=1:Nx+1
 for j=1:Ny+1
     if(f(j,i) >= 0)
       fplus(j,i)=f(j,i);
       fminus(j,i)=0;
     else
       fplus(j,i)=0;
       fminus(j,i)=-f(j,i);  
     end
 end
end
%Big matrixes
MX  = kron(Mx,speye(Ny+1));
MY  = kron(speye(Nx+1),My);

fintxy = sqrt(fplus(:))'*MX*MY*sqrt(fplus(:)) - ...
         sqrt(fminus(:))'*MX*MY*sqrt(fminus(:));
