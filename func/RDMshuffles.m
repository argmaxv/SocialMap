function ShfData=RDMshuffles(OrigData,n)

% With elements of OrigData, this function generate n suffled RDM while
% keeping its structure (mirrored imated)
% This function assume that all diag(OrigData) is 0
    
    xi=size(OrigData,1);
    mask_rdm=triu(ones(xi,xi),1);
    mask_mask=[1:1:xi*xi];
    Model_rdm_triu=OrigData(nonzeros(mask_mask.*mask_rdm(:)'));

    if nargin==1
        n=1;
    end
    m=0;
    while m<n
        clear X ShuffledData
        ShuffledData=zeros(xi,xi);
        X=randperm(numel(Model_rdm_triu));
        Y=0;
        for i=1:xi*xi
            if mask_rdm(i)==0
                ShuffledData(i)=0;
            else
                Y=Y+1;
                ShuffledData(i)=Model_rdm_triu(X(Y));
            end
        end
        clear ShuffledMat;
        ShuffledMat=triu(ShuffledData,1)'+ShuffledData;
        
        %ShuffledData=reshape(OrigData(X),size(OrigData));
        if OrigData==ShuffledMat
        else
            m=m+1;
            ShfData{m,1}=ShuffledMat;
        end
    end
    
end
