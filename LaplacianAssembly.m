% Generate Laplacian operator matrix B under moleduled slip boundary
% A+\lambda B = 0 here A is biharmonic operator and B is laplacian operator
% 2nd Order Accuracy !!

function [B] = LaplacianAssembly(L,H,Nx,Ny)

Dx = 2*L/(Nx-1);
Dy = 2*H/(Ny-1);

DX2i = 1/Dx/Dx;
DY2i = 1/Dy/Dy;

LaplacianMatrix = [0    ,          DY2i             ,   0   ;
                                DX2i,   -2*(DX2i+DY2i), DX2i;
                                0      ,        DY2i            ,     0  ];

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
% Here we start to assemble the Matrix B
N = (Nx-1)*(Ny-2);
MatB = zeros(N,N);

% j = 2 (y = -H + Dy) total 4 elements not 0 
for i = 1:1:Nx-1
    if i == 1
        rowIdx = PhiView21(Nx,i,2);
        MatB(rowIdx,1) = MatLaplacianAt(LaplacianMatrix,0,0);
        MatB(rowIdx,2) = MatLaplacianAt(LaplacianMatrix,1,0);
        MatB(rowIdx,Nx-1) = MatLaplacianAt(LaplacianMatrix,-1,0);
        MatB(rowIdx,Nx) = MatLaplacianAt(LaplacianMatrix,0,1);
    elseif i == Nx-1
        rowIdx = PhiView21(Nx,i,2);
        MatB(rowIdx,1) = MatLaplacianAt(LaplacianMatrix,1,0);
        MatB(rowIdx,Nx-2) = MatLaplacianAt(LaplacianMatrix,-1,0);
        MatB(rowIdx,Nx-1) = MatLaplacianAt(LaplacianMatrix,0,0);
        MatB(rowIdx,2*Nx-2) = MatLaplacianAt(LaplacianMatrix,0,1);
    else
        rowIdx = PhiView21(Nx,i,2);
        MatB(rowIdx,i-1) = MatLaplacianAt(LaplacianMatrix,-1,0);
        MatB(rowIdx,i) = MatLaplacianAt(LaplacianMatrix,0,0);
        MatB(rowIdx,i+1) = MatLaplacianAt(LaplacianMatrix,1,0);
        MatB(rowIdx,i+Nx-1) = MatLaplacianAt(LaplacianMatrix,0,1);
    end
end

% j = Ny-1 (y = H - Dy) total 4 elements not 0 
for i = 1:1:Nx-1
    if i == 1
        rowIdx = PhiView21(Nx,i,Ny-1);
        MatB(rowIdx,(Ny-4)*(Nx-1)+1) = MatLaplacianAt(LaplacianMatrix,0,-1);
        MatB(rowIdx,(Ny-3)*(Nx-1)+1) = MatLaplacianAt(LaplacianMatrix,0,0);
        MatB(rowIdx,(Ny-3)*(Nx-1)+2) = MatLaplacianAt(LaplacianMatrix,1,0);
        MatB(rowIdx,(Ny-3)*(Nx-1)+Nx-1) = MatLaplacianAt(LaplacianMatrix,-1,0);
    elseif i == Nx-1
        rowIdx = PhiView21(Nx,i,Ny-1);
        MatB(rowIdx,(Ny-4)*(Nx-1)+Nx-1) = MatLaplacianAt(LaplacianMatrix,0,-1);
        MatB(rowIdx,(Ny-3)*(Nx-1)+1) = MatLaplacianAt(LaplacianMatrix,1,0);
        MatB(rowIdx,(Ny-3)*(Nx-1)+Nx-2) = MatLaplacianAt(LaplacianMatrix,-1,0);
        MatB(rowIdx,(Ny-3)*(Nx-1)+Nx-1) = MatLaplacianAt(LaplacianMatrix,0,0);
    else
        rowIdx = PhiView21(Nx,i,Ny-1);
        MatB(rowIdx,(Ny-4)*(Nx-1)+i) = MatLaplacianAt(LaplacianMatrix,0,-1);
        MatB(rowIdx,(Ny-3)*(Nx-1)+i-1) = MatLaplacianAt(LaplacianMatrix,-1,0);
        MatB(rowIdx,(Ny-3)*(Nx-1)+i) = MatLaplacianAt(LaplacianMatrix,0,0);
        MatB(rowIdx,(Ny-3)*(Nx-1)+i+1) = MatLaplacianAt(LaplacianMatrix,1,0);
    end
end

% j >= 3 && j <= Ny-2  (total 5 elements not 0)
outFreq = floor(Ny/10);
for j = 3:1:Ny-2
    if mod(j,outFreq) == 0
        c = j/outFreq*10;
        disp(['Assemble Process: ' , num2str(c) , '%']);
    end
    for i = 1:1:Nx-1
        rowIdx = PhiView21(Nx,i,j);
        if i == 1
            MatB(rowIdx,(j-3)*(Nx-1)+1) = MatLaplacianAt(LaplacianMatrix,0,-1);
            MatB(rowIdx,(j-2)*(Nx-1)+1) = MatLaplacianAt(LaplacianMatrix,0,0);
            MatB(rowIdx,(j-2)*(Nx-1)+2) = MatLaplacianAt(LaplacianMatrix,1,0);
            MatB(rowIdx,(j-2)*(Nx-1)+Nx-1) = MatLaplacianAt(LaplacianMatrix,-1,0);
            MatB(rowIdx,(j-1)*(Nx-1)+1) = MatLaplacianAt(LaplacianMatrix,0,1);
        elseif i == Nx-1
            MatB(rowIdx,(j-3)*(Nx-1)+Nx-1) = MatLaplacianAt(LaplacianMatrix,0,-1);
            MatB(rowIdx,(j-2)*(Nx-1)+1) = MatLaplacianAt(LaplacianMatrix,1,0);
            MatB(rowIdx,(j-2)*(Nx-1)+Nx-2) = MatLaplacianAt(LaplacianMatrix,-1,0);
            MatB(rowIdx,(j-2)*(Nx-1)+Nx-1) = MatLaplacianAt(LaplacianMatrix,0,0);
            MatB(rowIdx,(j-1)*(Nx-1)+Nx-1) = MatLaplacianAt(LaplacianMatrix,0,1);
        else
            MatB(rowIdx,(j-3)*(Nx-1)+i) = MatLaplacianAt(LaplacianMatrix,0,-1);
            MatB(rowIdx,(j-2)*(Nx-1)+i-1) = MatLaplacianAt(LaplacianMatrix,-1,0);
            MatB(rowIdx,(j-2)*(Nx-1)+i) = MatLaplacianAt(LaplacianMatrix,0,0);
            MatB(rowIdx,(j-2)*(Nx-1)+i+1) = MatLaplacianAt(LaplacianMatrix,1,0);
            MatB(rowIdx,(j-1)*(Nx-1)+i) = MatLaplacianAt(LaplacianMatrix,0,1);
        end
    end
end

% Laplacian Matrix B has been assembled
% Check its symmetry property
Error = norm(MatB-MatB',1);
if Error == 0
    fprintf('Laplacian Matrix B has been assembled with error = %f \n', Error);
end

B = MatB;
clear MatB;

end