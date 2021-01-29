function [E_rdm, Ctx_rdm, Gr_rdm, D_rdm, I_rdm, E_rdm_perm, D_rdm_perm, I_rdm_perm]=EucRDM_factorize_ss(template, Nperm)
    
    rng('default');
    sn_size=16; %social network size (n individuals)
    q=[1:4]; %4x4 social hierarchy structure

    for jdx=1:size(template,1)
        contexts(jdx,1)=str2double(template{jdx}(2)); %context
        faces(jdx,1)=str2double(template{jdx}(4:5)); % face id
    end
    CtxMask=squareform(pdist(contexts, 'euclidean')); Ctx_rdm=CtxMask>0;
    [~, Grp, ~]=faceinfo(faces);
    GrMask=squareform(pdist(Grp', 'euclidean')); Gr_rdm=GrMask>0;

    [Xo,Yo]=meshgrid(q,q);
    X=Xo(faces);
    Y=Yo(faces);    
    Rel=[X(contexts==1); Y(contexts==2)];
    Irr =[Y(contexts==1); X(contexts==2)];

    I=[1:length(X)];
    [i,j]=meshgrid(I,I);

    E_rdm=sqrt((X(i)-X(j)).^2+(Y(i)-Y(j)).^2);
    D_rdm=abs(Rel(i)-Rel(j));
    I_rdm=abs(Irr(i)-Irr(j));
    
    if nargin>1 %Permutation
        Xperm=randi([1,4],sn_size,Nperm);
        Yperm=randi([1,4],sn_size,Nperm);        
        for n=1:Nperm
            clear X0 Y0
            X0_coord=Xperm(:,n);
            Y0_coord=Yperm(:,n);
            X0=X0_coord(faces);
            Y0=Y0_coord(faces);
            D0=[X0(contexts==1); Y0(contexts==2)];
            I0=[Y0(contexts==1); X0(contexts==2)];
            I=[1:length(X)];
            [i,j]=meshgrid(I,I);
            E_rdm_perm{n,1}=sqrt((X0(i)-X0(j)).^2+(Y0(i)-Y0(j)).^2);
            D_rdm_perm{n,1}=abs(D0(i)-D0(j));
            I_rdm_perm{n,1}=abs(I0(i)-I0(j));
        end %for Nperm
    end %if permutation

end
    