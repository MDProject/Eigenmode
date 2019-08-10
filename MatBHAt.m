% return the element of matrix BiharMatrix at position (i,j)
% (i,j) is the center-referenced index of {-2,-1,0,1,2}*{-2,-1,0,1,2}
% undefined elements is set to default value 0
function Cij = MatBHAt(C,i,j)
    ci0 = 3;
    cj0 = 3;
    Cij = C(ci0+i,cj0+j);
end
