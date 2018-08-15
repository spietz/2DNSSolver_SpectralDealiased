function [fintxy] = integrate2d(f, Nx, Ny, Mx, My)
% remember we always integrate |integrant|² using discrete matrix

% integrate x dir (by rows)
fintx = zeros(Ny+1,1);
for i = 1:Ny+1
    fintx(i) = sqrt(f(i,:))*(Mx)*sqrt(f(i,:))';
end
% integrate y dir (the column)
fintxy = sqrt(fintx)'*(My)*sqrt(fintx);
