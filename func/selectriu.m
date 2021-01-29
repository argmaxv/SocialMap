
function [mtxur]=selectriu(mtx) 
% Averaging the upper-right and lower-left sides of the matrix and returns the values of one side of matix

    [i,j]=size(mtx);
    if i~=j
        error('Check the matrix size.')
    end
    mtx_symmetry=(mtx+mtx')/2;
    temp=triu(ones(i,j));
    temp(temp==0)=NaN;
    mtxur0=mtx_symmetry.*temp;
    mtxur=rmmissing(mtxur0(:));
    
end
