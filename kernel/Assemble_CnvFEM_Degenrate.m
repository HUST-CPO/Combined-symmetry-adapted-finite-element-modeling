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
    epsilonr=physic.epsilonr(domain);
    mu=physic.mur(domain);
    
    %submatrix
    St=zeros(3,3);Tte=zeros(3,3);
    Sz=zeros(3,3);Tz=zeros(3,3);
    G=zeros(3,3);Ttu=zeros(3,3);
    for i=1:3
       for j=1:3
          for k=1:IntegralOrder
              St(i,j)=St(i,j)+pt(k)*DetJac/mu*dot(curlEt(i,:,k),curlEt(j,:,k));
              Tte(i,j)=Tte(i,j)+pt(k)*DetJac*k0*k0*epsilonr*dot(Et(i,:,k),Et(j,:,k));
              Sz(i,j)=Sz(i,j)+pt(k)*DetJac/mu*dot(gradEz(i,:,k),gradEz(j,:,k));
              Tz(i,j)=Tz(i,j)+pt(k)*DetJac*k0*k0*epsilonr*dot(Ez(i,:,k),Ez(j,:,k));
              G(i,j)=G(i,j)+pt(k)*DetJac/mu*dot(Et(i,:,k),gradEz(j,:,k));
              Ttu(i,j)=Ttu(i,j)+pt(k)*DetJac/mu*dot(Et(i,:,k),Et(j,:,k));
          end
       end
    end
    Gt=G';
    
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
            Bv(NumNonB)=G(i,j);

            NumNonB=NumNonB+1;
            Bi(NumNonB)=MappingIndexSi;
            Bj(NumNonB)=MappingIndexVj;
            Bv(NumNonB)=Gt(i,j);
        end
    end
end

A = sparse(Ai,Aj,Av);
B = sparse(Bi,Bj,Bv);
A=sparse([A zeros(mesh.Dof);zeros(mesh.Dof) A]);
B=sparse([B zeros(mesh.Dof);zeros(mesh.Dof) B]);

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

Index=sort(unique([mesh.BBCIndex+mesh.Dof;mesh.PECIndex;mesh.PECIndex2+mesh.Dof]),'descend');
for i=1:length(Index)
    tempN=Index(i);
    P(:,tempN)=[];
end

solver.P=P;

solver.A=P'*A*P;
solver.B=P'*B*P;

end