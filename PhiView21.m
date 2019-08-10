% transform the 2D index (i,j) of Phi into 1D according to the natural
% sequence
%           *****(y = H, j = Ny)*********************
%           **(i: 1~Nx-1, j = Ny-1)******************---- 
%           ****************************************  |
%                                                                                |
%                                                                          Phi inner points   
%                                                                                |
%           ****************************************  |
%           **(i: 1~Nx-1, j = 2)**********************---
%           *****(y = -H, j = 1)**********************

function I = PhiView21(Nx,i,j)
    I = (j-2)*(Nx-1) + i;
end
