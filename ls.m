% return the slip length at each grid point
% idx:  grid point index according to inertial phi points natural sequence 
% idx=0 denotes the position at (-L,-H)
function slipLength = ls(idx, DX, L)
    l0 = 1;
    kx = pi/L;
    slipLength = l0 + sin(kx*(-L+(idx-1)*DX)); 
end
