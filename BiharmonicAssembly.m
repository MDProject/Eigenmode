% Generate Biharmonic differential matrix A under moleduled slip boundary
% A+\lambda B = 0 here A is biharmonic operator and B is laplacian operator

function [A] = BiharmonicAssembly(L,H,Nx,Ny)
% Nx is total grid points along x direction(start&end points included)
% Ny is total grid points along y direction(start&end points included)
% L is half length along x direction
% H is half height along y direction
Dx = 2*L/(Nx-1);
Dy = 2*H/(Ny-1);

DX2i = 1/Dx/Dx;
DY2i = 1/Dy/Dy;
DX4i = DX2i^2;
DY4i = DY2i^2;
DX2DY2i = DX2i*DY2i;

% 13 points 4th-order biharmonic difference scheme
BiharMatrix = [ 0        ,          0                      ,         DY4i                              ,         0                     ,   0    ;
                           0        ,   2*DX2DY2i            ,-4*DY4i-4*DX2i*DY2i           ,     2*DX2DY2i         ,   0    ;
                           DX4i  ,-4*DX4i-4*DX2DY2i,6*DY4i+6*DX4i+8*DX2DY2i,-4*DX4i-4*DX2DY2i,DX4i;
                           0        ,   2*DX2DY2i            ,-4*DY4i-4*DX2DY2i              ,     2*DX2DY2i         ,  0    ;
                           0        ,       0                         ,         DY4i                              ,         0                     ,  0    ];

% Ghost-points ratio, Ri in the note
% Phi ~ (1~Nx-1)*(2~Ny-1) , scalar potential points to be solved; ghost
% points index domain [1~Nx-1]
Ratio = zeros(1,Nx-1);
for i = 1:1:Nx-1
    Ratio(i) = (Dy - 2*ls(i,Dx,L))/(Dy + 2*ls(i,Dx,L));
end

% Scalar potential Phi ~ (1~Nx-1)*(2~Ny-1), which is aligned according to
% the note is a column vector of size (1~Nx-1)*(2~Ny-1) since Phi is fixed
% to be 0 at j = 1 and Ny
% Here we start to assemble the Matrix A
N = (Nx-1)*(Ny-2);
MatA = zeros(N,N);
% j = 2 (y = -H) total 9 elements not 0 
for i = 1:1:Nx-1
    if i ==1
        rowIdx = PhiView21(Nx,i,2);
        MatA(rowIdx,1) = MatBHAt(BiharMatrix,0,-2)*Ratio(1) + MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,2) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,3) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,Nx-2) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,Nx-1) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,Nx) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,Nx+1) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,2*Nx-2) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,2*Nx-1) = MatBHAt(BiharMatrix,0,2);
    elseif i == 2
        rowIdx = PhiView21(Nx,2,2);
        MatA(rowIdx,1) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,2) = MatBHAt(BiharMatrix,0,-2)*Ratio(2) + MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,3) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,4) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,Nx-1) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,Nx) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,Nx+1) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,Nx+2) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,2*Nx) = MatBHAt(BiharMatrix,0,2);
    elseif i == Nx-2
        rowIdx = PhiView21(Nx,i,2);
        MatA(rowIdx,1) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,Nx-4) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,Nx-3) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,Nx-2) = MatBHAt(BiharMatrix,0,-2)*Ratio(Nx-2) + MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,Nx-1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,2*Nx-4) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,2*Nx-3) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,2*Nx-2) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,3*Nx-4) = MatBHAt(BiharMatrix,0,2);
    elseif i == Nx-1
        rowIdx = PhiView21(Nx,i,2);
        MatA(rowIdx,1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,2) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,Nx-3) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,Nx-1) = MatBHAt(BiharMatrix,0,-2)*Ratio(Nx-1) + MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,Nx-2) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,Nx) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,2*Nx-3) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,2*Nx-2) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,3*Nx-3) = MatBHAt(BiharMatrix,0,2);
    else
        rowIdx = PhiView21(Nx,i,2);
        MatA(rowIdx,i-2) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,i-1) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,i) = MatBHAt(BiharMatrix,0,-2)*Ratio(i) + MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,i+1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,i+2) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,Nx+i-2) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,Nx+i-1) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,Nx+i) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,2*Nx+i-2) = MatBHAt(BiharMatrix,0,2);
    end
end
% j = 3 total 12 elements not 0 
for i = 1:1:Nx-1
    if i == 1
        rowIdx = PhiView21(Nx,i,3);
        MatA(rowIdx,1) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,2) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,Nx-1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,Nx) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,Nx+1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,Nx+2) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,2*Nx-3) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,2*Nx-2) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,2*Nx-1) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,2*Nx) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,3*Nx-3) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,3*Nx-2) = MatBHAt(BiharMatrix,0,2);
    elseif i ==2
        rowIdx = PhiView21(Nx,i,3);
        MatA(rowIdx,1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,2) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,3) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,Nx) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,Nx+1) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,Nx+2) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,Nx+3) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,2*Nx-2) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,2*Nx-1) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,2*Nx) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,2*Nx+1) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,3*Nx-1) = MatBHAt(BiharMatrix,0,2);
    elseif i == Nx-2
        rowIdx = PhiView21(Nx,i,3);
        MatA(rowIdx,Nx-3) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,Nx-2) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,Nx-1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,Nx) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,2*Nx-5) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,2*Nx-4) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,2*Nx-3) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,2*Nx-2) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,3*Nx-5) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,3*Nx-4) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,3*Nx-3) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,4*Nx-5) = MatBHAt(BiharMatrix,0,2);
    elseif i == Nx-1
        rowIdx = PhiView21(Nx,i,3);
        MatA(rowIdx,1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,Nx-2) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,Nx-1) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,Nx) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,Nx+1) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,2*Nx-4) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,2*Nx-3) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,2*Nx-2) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,2*Nx-1) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,3*Nx-4) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,3*Nx-3) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,4*Nx-4) = MatBHAt(BiharMatrix,0,2);
    else
        rowIdx = PhiView21(Nx,i,3);
        MatA(rowIdx,i-1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,i) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,i+1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,Nx+i-3) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,Nx+i-2) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,Nx+i-1) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,Nx+i) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,Nx+i+1) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,2*Nx+i-3) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,2*Nx+i-2) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,2*Nx+i-1) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,3*Nx+i-3) = MatBHAt(BiharMatrix,0,2);
    end
end
% j = Ny-2 total 12 elements not 0 
for i = 1:1:Nx-1
    if i == 1
        rowIdx = PhiView21(Nx,i,Ny-2);
        MatA(rowIdx,(Ny-6)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-5)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+2) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+2) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+3) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+2) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,(Ny-2)*(Nx-1)) = MatBHAt(BiharMatrix,-1,1);
    elseif i == 2
        rowIdx = PhiView21(Nx,i,Ny-2);
        MatA(rowIdx,(Ny-6)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-5)*(Nx-1)+1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+3) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+1) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+3) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+4) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+3) = MatBHAt(BiharMatrix,1,1);
    elseif i == Nx-2
        rowIdx = PhiView21(Nx,i,Ny-2);
        MatA(rowIdx,(Ny-6)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-5)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+1) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-4) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,1,1);
    elseif i == Nx-1
        rowIdx = PhiView21(Nx,i,Ny-2);
        MatA(rowIdx,(Ny-6)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-5)*(Nx-1)+1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+2) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,1,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,1);
    else
        rowIdx = PhiView21(Nx,i,Ny-2);
        MatA(rowIdx,(Ny-6)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-5)*(Nx-1)+i-1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+i+1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+i-2) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+i-1) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+i+1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-4)*(Nx-1)+i+2) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+i-1) = MatBHAt(BiharMatrix,-1,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+i+1) = MatBHAt(BiharMatrix,1,1);
    end
end
% j = Ny-1 total 9 elements not 0 in each row
for i = 1:1:Nx-1
    if i ==1
        rowIdx = PhiView21(Nx,i,Ny-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-4)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+2) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,0) + MatBHAt(BiharMatrix,0,2)*Ratio(1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+2) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+3) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,-1,0);
    elseif i == 2
        rowIdx = PhiView21(Nx,i,Ny-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-4)*(Nx-1)+1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+3) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,0) + MatBHAt(BiharMatrix,0,2)*Ratio(2);
        MatA(rowIdx,(Ny-3)*(Nx-1)+3) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+4) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,-2,0);
    elseif i == Nx-2
        rowIdx = PhiView21(Nx,i,Ny-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-4) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,0) + MatBHAt(BiharMatrix,0,2)*Ratio(Nx-2);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,1,0);
    elseif i == Nx-1
        rowIdx = PhiView21(Nx,i,Ny-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-4)*(Nx-1)+1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+2) = MatBHAt(BiharMatrix,2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,0) + MatBHAt(BiharMatrix,0,2)*Ratio(Nx-1);
    else
        rowIdx = PhiView21(Nx,i,Ny-1);
        MatA(rowIdx,(Ny-5)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,-2);
        MatA(rowIdx,(Ny-4)*(Nx-1)+i-1) = MatBHAt(BiharMatrix,-1,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,-1);
        MatA(rowIdx,(Ny-4)*(Nx-1)+i+1) = MatBHAt(BiharMatrix,1,-1);
        MatA(rowIdx,(Ny-3)*(Nx-1)+i-2) = MatBHAt(BiharMatrix,-2,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+i-1) = MatBHAt(BiharMatrix,-1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,0) + MatBHAt(BiharMatrix,0,2)*Ratio(i);
        MatA(rowIdx,(Ny-3)*(Nx-1)+i+1) = MatBHAt(BiharMatrix,1,0);
        MatA(rowIdx,(Ny-3)*(Nx-1)+i+2) = MatBHAt(BiharMatrix,2,0);
    end
end
% j >=4 && j<=Ny-3 inertial layers, all 13 points are used
outFreq = floor(Ny/10);
for j = 4:1:Ny-3
    if mod(j,outFreq) == 0
        c = j/outFreq*10;
        disp(['Assemble Process: ' , num2str(c) , '%']);
    end
    for i = 1:1:Nx-1
        rowIdx = PhiView21(Nx,i,j);
        if i == 1
            MatA(rowIdx,(j-4)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,-2);
            MatA(rowIdx,(j-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+2) = MatBHAt(BiharMatrix,1,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,-1,-1);
            MatA(rowIdx,(j-2)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,0);
            MatA(rowIdx,(j-2)*(Nx-1)+2) = MatBHAt(BiharMatrix,1,0);
            MatA(rowIdx,(j-2)*(Nx-1)+3) = MatBHAt(BiharMatrix,2,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-2,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,-1,0);
            MatA(rowIdx,(j-1)*(Nx-1)+1) = MatBHAt(BiharMatrix,0,1);
            MatA(rowIdx,(j-1)*(Nx-1)+2) = MatBHAt(BiharMatrix,1,1);
            MatA(rowIdx,(j-1)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,-1,1);
            MatA(rowIdx, j*(Nx-1)+1) = MatBHAt(BiharMatrix,0,2);
        elseif i == 2
            MatA(rowIdx,(j-4)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,-2);
            MatA(rowIdx,(j-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,-1,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+3) = MatBHAt(BiharMatrix,1,-1);
            MatA(rowIdx,(j-2)*(Nx-1)+1) = MatBHAt(BiharMatrix,-1,0);
            MatA(rowIdx,(j-2)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,0);
            MatA(rowIdx,(j-2)*(Nx-1)+3) = MatBHAt(BiharMatrix,1,0);
            MatA(rowIdx,(j-2)*(Nx-1)+4) = MatBHAt(BiharMatrix,2,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,-2,0);
            MatA(rowIdx,(j-1)*(Nx-1)+1) = MatBHAt(BiharMatrix,-1,1);
            MatA(rowIdx,(j-1)*(Nx-1)+2) = MatBHAt(BiharMatrix,0,1);
            MatA(rowIdx,(j-1)*(Nx-1)+3) = MatBHAt(BiharMatrix,1,1);
            MatA(rowIdx, j*(Nx-1)+2) = MatBHAt(BiharMatrix,0,2);
        elseif i == Nx-2
            MatA(rowIdx,(j-4)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,-2);
            MatA(rowIdx,(j-3)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-1,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,1,-1);
            MatA(rowIdx,(j-2)*(Nx-1)+1) = MatBHAt(BiharMatrix,2,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-4) = MatBHAt(BiharMatrix,-2,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-1,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,1,0);
            MatA(rowIdx,(j-1)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-1,1);
            MatA(rowIdx,(j-1)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,1);
            MatA(rowIdx,(j-1)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,1,1);
            MatA(rowIdx, j*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,0,2);
        elseif i == Nx-1
            MatA(rowIdx,(j-4)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,-2);
            MatA(rowIdx,(j-3)*(Nx-1)+1) = MatBHAt(BiharMatrix,1,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-1,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,-1);
            MatA(rowIdx,(j-2)*(Nx-1)+1) = MatBHAt(BiharMatrix,1,0);
            MatA(rowIdx,(j-2)*(Nx-1)+2) = MatBHAt(BiharMatrix,2,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-3) = MatBHAt(BiharMatrix,-2,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-1,0);
            MatA(rowIdx,(j-2)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,0);
            MatA(rowIdx,(j-1)*(Nx-1)+1) = MatBHAt(BiharMatrix,1,1);
            MatA(rowIdx,(j-1)*(Nx-1)+Nx-2) = MatBHAt(BiharMatrix,-1,1);
            MatA(rowIdx,(j-1)*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,1);
            MatA(rowIdx, j*(Nx-1)+Nx-1) = MatBHAt(BiharMatrix,0,2);
        else
            MatA(rowIdx,(j-4)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,-2);
            MatA(rowIdx,(j-3)*(Nx-1)+i-1) = MatBHAt(BiharMatrix,-1,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,-1);
            MatA(rowIdx,(j-3)*(Nx-1)+i+1) = MatBHAt(BiharMatrix,1,-1);
            MatA(rowIdx,(j-2)*(Nx-1)+i-2) = MatBHAt(BiharMatrix,-2,0);
            MatA(rowIdx,(j-2)*(Nx-1)+i-1) = MatBHAt(BiharMatrix,-1,0);
            MatA(rowIdx,(j-2)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,0);
            MatA(rowIdx,(j-2)*(Nx-1)+i+1) = MatBHAt(BiharMatrix,1,0);
            MatA(rowIdx,(j-2)*(Nx-1)+i+2) = MatBHAt(BiharMatrix,2,0);
            MatA(rowIdx,(j-1)*(Nx-1)+i-1) = MatBHAt(BiharMatrix,-1,1);
            MatA(rowIdx,(j-1)*(Nx-1)+i) = MatBHAt(BiharMatrix,0,1);
            MatA(rowIdx,(j-1)*(Nx-1)+i+1) = MatBHAt(BiharMatrix,1,1);
            MatA(rowIdx, j*(Nx-1)+i) = MatBHAt(BiharMatrix,0,2);
        end
    end
end
% Biharmonic Matrix A has been assembled
% Check its symmetry property
Error = norm(MatA-MatA',1);
if Error == 0
    fprintf('Biharmonic Matrix A has been assembled with error = %f \n', Error);
end

A = MatA;
clear MatA;

end