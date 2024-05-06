function [solver,mesh]=Assemble_CnvFEM_Degenrate(mesh,physic)

k0=physic.k0;
%guass points
[xt,yt,pt,IntegralOrder]=GetGuassPoints(2);

%init
mesh.Dof=mesh.nbrVertex+mesh.nbrEdges;
NbrNon=mesh.nbrTri*9;
NumNonA=0;
Ai=zeros(NbrNon,1);
Aj=zeros(NbrNon,1);
Av=zeros(NbrNon,1);
NumNonB=0;
Bi=zeros(NbrNon*4,1);
Bj=zeros(NbrNon*4,1);
Bv=zeros(NbrNon*4,1);

%loop tri
for n=1:mesh.nbrTri
    %coordinate of tri's vertex
    x=zeros(3,1);y=zeros(3,1);l=zeros(3,1);
    x(1)=mesh.vertex(mesh.tri(n,1),1);y(1)=mesh.vertex(mesh.tri(n,1),2);
    x(2)=mesh.vertex(mesh.tri(n,2),1);y(2)=mesh.vertex(mesh.tri(n,2),2);
    x(3)=mesh.vertex(mesh.tri(n,3),1);y(3)=mesh.vertex(mesh.tri(n,3),2);
    l(1)=sqrt((x(1)-x(2))*(x(1)-x(2))+(y(1)-y(2))*(y(1)-y(2)));
    l(2)=sqrt((x(1)-x(3))*(x(1)-x(3))+(y(1)-y(3))*(y(1)-y(3)));
    l(3)=sqrt((x(3)-x(2))*(x(3)-x(2))+(y(3)-y(2))*(y(3)-y(2)));

    %jac matrix
    Jac=zeros(3,3);
    Jac(1,1)=x(2)-x(1);Jac(1,2)=y(2)-y(1);
    Jac(2,1)=x(3)-x(1);Jac(2,2)=y(3)-y(1);Jac(3,3)=1;
    InvJac=inv(Jac);
    DetJac=abs(det(Jac));
    DetJacx=det(Jac);
    TJac=Jac'/DetJacx;

    %bf
    Et=zeros(3,3,IntegralOrder);curlEt=zeros(3,3,IntegralOrder);
    Ez=zeros(3,3,IntegralOrder);gradEz=zeros(3,3,IntegralOrder);
    temp=zeros(3,1);
    for i=1:IntegralOrder
        for j=1:3
            [temp(1),temp(2)]=BF_Et(j,xt(i),yt(i));temp(3)=0;
            Et(j,:,i)=InvJac*temp*l(j);
            temp(3)=BF_curlEt(j,xt(i),yt(i));temp(1)=0;temp(2)=0;
            curlEt(j,:,i)=TJac*temp*l(j);
            temp(3)=BF_Ez(j,xt(i),yt(i));temp(1)=0;temp(2)=0;
            Ez(j,:,i)=temp;
            [temp(1),temp(2)]=BF_gradEz(j,xt(i),yt(i));temp(3)=0;
            gradEz(j,:,i)=InvJac*temp;
        end
    end

    %material
    domain=mesh.triID(n);
    epsr=eye(3)*physic.epsilonr(domain);
    mur=eye(3)*physic.mur(domain);
    %PML 
    R0=100;
    if ismember(domain,physic.PML)
        r0=physic.PMLData(1);
        dr=physic.PMLData(2);

        averX=(x(1)+x(2)+x(3))/3;
        averY=(y(1)+y(2)+y(3))/3;
        [~,rho]=cart2pol(averX,averY);
        sigma=((rho-r0)/dr)^2*R0/dr/k0;
        s1=1-1i*dr/2/rho*sigma;
        s2=1-1i*sigma;
        aa=s1/s2;bb=s2/s1;cc=s1*s2;
        Lambda=[(aa*averX*averX+bb*averY*averY)/rho/rho,((aa-bb)*averX*averY)/rho/rho,0;
                ((aa-bb)*averX*averY)/rho/rho,(bb*averX*averX+aa*averY*averY)/rho/rho,0;
                0,0,cc];
        epsr=epsr*Lambda;
        mur=mur*Lambda;
    end
    murt=zeros(3,3);
    mur=inv(mur);
    epsz=epsr(3,3);
    epst=epsr;
    epst(3,3)=0;
    murz=mur(3,3);
    murt(1,1)=mur(2,2);murt(1,2)=-mur(2,1);
    murt(2,1)=-mur(1,2);murt(2,2)=mur(1,1);
    
    %submatrix
    St=zeros(3,3);Tte=zeros(3,3);
    Sz=zeros(3,3);Tz=zeros(3,3);
    G=zeros(3,3);Ttu=zeros(3,3);
    Gt=zeros(3,3);
    for i=1:3
       for j=1:3
          for k=1:IntegralOrder
              St(i,j)=St(i,j)+pt(k)*DetJac*murz*dot(curlEt(i,:,k),curlEt(j,:,k));
              Tte(i,j)=Tte(i,j)+pt(k)*DetJac*k0*k0*sum(Et(i,:,k)'.*(epst*Et(j,:,k)'));
              Sz(i,j)=Sz(i,j)+pt(k)*DetJac*sum((murt*gradEz(j,:,k)').*gradEz(i,:,k)');
              Tz(i,j)=Tz(i,j)+pt(k)*DetJac*k0*k0*epsz*dot(Ez(i,:,k),Ez(j,:,k));
              G(i,j)=G(i,j)+pt(k)*DetJac*sum((murt*Et(j,:,k)').*gradEz(i,:,k)');
              Gt(i,j)=Gt(i,j)+pt(k)*DetJac*sum((murt*gradEz(j,:,k)').*Et(i,:,k)');
              Ttu(i,j)=Ttu(i,j)+pt(k)*DetJac*sum((murt*Et(j,:,k)').*Et(i,:,k)');
          end
       end
    end
    
    % 排入三元组
    for i=1:3
        for j=1:3
            MappingIndexVi=mesh.edgesOfTri(n,i)+mesh.nbrVertex;
            MappingIndexVj=mesh.edgesOfTri(n,j)+mesh.nbrVertex;
            MappingIndexSi=mesh.tri(n,i);
            MappingIndexSj=mesh.tri(n,j);
            
            NumNonA=NumNonA+1;
            Ai(NumNonA)=MappingIndexVi;
            Aj(NumNonA)=MappingIndexVj;
            Av(NumNonA)=St(i,j)-Tte(i,j);

            NumNonB=NumNonB+1;
            Bi(NumNonB)=MappingIndexVi;
            Bj(NumNonB)=MappingIndexVj;
            Bv(NumNonB)=Ttu(i,j);

            NumNonB=NumNonB+1;
            Bi(NumNonB)=MappingIndexSi;
            Bj(NumNonB)=MappingIndexSj;
            Bv(NumNonB)=Sz(i,j)-Tz(i,j);

            NumNonB=NumNonB+1;
            Bi(NumNonB)=MappingIndexVi;
            Bj(NumNonB)=MappingIndexSj;
            Bv(NumNonB)=Gt(i,j);

            NumNonB=NumNonB+1;
            Bi(NumNonB)=MappingIndexSi;
            Bj(NumNonB)=MappingIndexVj;
            Bv(NumNonB)=G(i,j);
        end
    end
end

A = sparse(Ai,Aj,Av);
B = sparse(Bi,Bj,Bv);
A=sparse([A sparse(mesh.Dof,mesh.Dof);sparse(mesh.Dof,mesh.Dof) A]);
B=sparse([B sparse(mesh.Dof,mesh.Dof);sparse(mesh.Dof,mesh.Dof) B]);

P=speye(mesh.Dof*2);
%treat boundary condition
%Cnv boundary condition
phi=physic.phi*physic.m;
phi=tan(phi/2);
NbrBBC=length(mesh.BBCIndex);
for i=1:NbrBBC
    tempN=mesh.BBCIndex(i);
    P(mesh.Dof+tempN,tempN)=phi;
end

Index=sort(unique([1;mesh.BBCIndex+mesh.Dof;mesh.PECIndex;mesh.PECIndex+mesh.Dof;mesh.PECIndex2+mesh.Dof]),'descend');
for i=1:length(Index)
    tempN=Index(i);
    P(:,tempN)=[];
end

solver.P=P;

solver.A=P'*A*P;
solver.B=P'*B*P;

end