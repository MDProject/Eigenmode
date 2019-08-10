% Transform the inner points Phi into 2D array
% automatically pick up the points at y = +- H and set them to be 0 as well
% as fill up the points at x = +L

%   *********************************


function Phi = PhiReshape(P,Nx,Ny)

Phi = zeros(Ny,Nx);
for i = Ny:-1:1
    if i >=2 && i <= Ny-1
        Phi(i,1:Nx-1) = P((Ny-i-1)*(Nx-1)+1 : (Ny-i-1)*(Nx-1)+Nx-1);
    end
end
% link the periodic condition along x
Phi(:,Nx) = Phi(:,1);

end

% Test Sample
% P = linspace(1,20,20);
% P2D = PhiReshape(P,6,6);
