% Calculate the gradient of Phi (velocity field)
% [Phi(i-1,j) - Phi(i+1,j)]/2Dy = grad_y(vx) of Phi  and  [Phi(i,j+1) - Phi(i,j-1)]/2Dx = grad_x(-vy) pf Phi
% vx:   [Phi(i-1,j) - Phi(i+1,j)]/2Dy   vy:        [Phi(i,j-1) - Phi(i,j+1)]/2Dx
% shape of Phi: Ny*Nx

function [Vx,Vy] = GradientPhi(Phi,Nx,Ny,L,H)
Dx = 2*L/(Nx-1);
Dy = 2*H/(Ny-1);

% Used for constructing the ghost points
Ratio = zeros(1,Nx);
for i = 1:1:Nx
    Ratio(i) = (Dy - 2*ls(i,Dx,L))/(Dy + 2*ls(i,Dx,L));
end

PhiGhost0J = Ratio.*Phi(2,:);
PhiGhostNyplus1 = Ratio.*Phi(Ny-1,:); 
% Vx component
Vx = zeros(Ny,Nx);
for i = 2:1:Ny-1
    for j = 1:1:Nx
        Vx(i,j) = (Phi(i-1,j) - Phi(i+1,j))/2/Dy;
    end
end
for j = 1:1:Nx
    Vx(1,j) = (PhiGhost0J(j) - Phi(2,j))/2/Dy;
    Vx(Ny,j) = (Phi(Ny-1,j) - PhiGhostNyplus1(j))/2/Dy;
end


% Vy component
Vy = zeros(Ny,Nx);
for j = 2:1:Nx-1
    for i = 1:1:Ny
        Vy(i,j) = (Phi(i,j-1) - Phi(i,j+1))/2/Dx;
    end
end
for i = 1:1:Ny
    Vy(i,1) = (Phi(i,Nx-1) - Phi(i,2))/2/Dx;
    Vy(i,Nx) = (Phi(i,Nx-1) - Phi(i,2))/2/Dx;
end

end