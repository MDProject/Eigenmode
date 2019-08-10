% return the element of matrix BiharMatrix at position (i,j)
% (i,j) is the center-referenced index of {-1,0,1}*{-1,0,1}
% undefined elements is set to default value 0

function Cij = MatLaplacianAt(C,i,j)
    ci0 = 2;
    cj0 = 2;
    Cij = C(ci0+i,cj0+j);
end