function [solver,mesh]=Assemble_CnFEM(mesh,physic)

k0=physic.k0;
Lambda=2*pi/physic.N;
%guass points
[xt,yt,pt,IntegralOrder]=GetGuassPoints(2);

%PEC
boundaryFlag=mesh.PECID;
index=find(mesh.BedgeID==boundaryFlag(1));
for i=2:length(boundaryFlag)
    index=[index;find(mesh.BedgeID==boundaryFlag(i))];
end
index=sort(index);
Boundary=mesh.Bedge(index,:);
Boundary=sort(Boundary,2);
[PECNode,PECEdge]=GetPECIndex(Boundary,mesh.edges);
%rotation periodic boundary
NbrCn=length(mesh.SrcBoundary);
[CnzSrcTempIndex,CntSrcTempIndex]=GetPECIndex(mesh.SrcBoundary,mesh.edges);%src
[CnzDstIndex,CntDstIndex]=GetPECIndex(mesh.DstBoundary,mesh.edges);%dst
CnzSrcIndex=zeros(NbrCn+1,1);CntSrcIndex=zeros(NbrCn,1);
CnSign=zeros(NbrCn,1);
for i=1:NbrCn+1
    tempx1=mesh.vertex(CnzDstIndex(i),1);
    tempy1=mesh.vertex(CnzDstIndex(i),2);
    tempr1=sqrt(tempx1^2+tempy1^2);
    for j=1:NbrCn+2-i
        tempx2=mesh.vertex(CnzSrcTempIndex(j),1);
        tempy2=mesh.vertex(CnzSrcTempIndex(j),2);
        tempr2=sqrt(tempx2^2+tempy2^2);
        if abs(tempr2-tempr1)<4e-10
            CnzSrcIndex(i)=CnzSrcTempIndex(j);
            CnzSrcTempIndex(j)=[];
            break;
        end
    end
end
for i=1:NbrCn 
    tempnode1=mesh.edges(CntDstIndex(i),1);
    tempnode2=mesh.edges(CntDstIndex(i),2);
    tempnode1x=CnzSrcIndex(CnzDstIndex==tempnode1);
    tempnode2x=CnzSrcIndex(CnzDstIndex==tempnode2);
    if tempnode1x>tempnode2x 
        temp=tempnode1x;
        tempnode1x=tempnode2x;
        tempnode2x=temp;
    end
    tempx1=mesh.vertex(tempnode1,1);tempx2=mesh.vertex(tempnode2,1);
    tempy1=mesh.vertex(tempnode1,2);tempy2=mesh.vertex(tempnode2,2);
    tempr1=sqrt(tempx1^2+tempy1^2);tempr2=sqrt(tempx2^2+tempy2^2);
    if tempr1>tempr2
       sign1=-1;
    else
        sign1=1;
    end
    for j=1:NbrCn
        tempnode1=mesh.edges(CntSrcTempIndex(j),1);
        tempnode2=mesh.edges(CntSrcTempIndex(j),2);
        if (tempnode1==tempnode1x)&&(tempnode2==tempnode2x)
            CntSrcIndex(i)=CntSrcTempIndex(j);
            tempx1=mesh.vertex(tempnode1,1);tempx2=mesh.vertex(tempnode2,1);
            tempy1=mesh.vertex(tempnode1,2);tempy2=mesh.vertex(tempnode2,2);
            tempr1=sqrt(tempx1^2+tempy1^2);tempr2=sqrt(tempx2^2+tempy2^2);
            if tempr1>tempr2
                sign2=-1;
            else
                sign2=1;
            end
            CnSign(i)=sign1*sign2;
            break;
        end
    end
end
CnzDstIndex(CnzDstIndex==mesh.OrigIndex)=[];
CnzSrcIndex(CnzSrcIndex==mesh.OrigIndex)=[];

%init
mesh.Dof=mesh.nbrEdges+mesh.nbrVertex;
NbrNon=mesh.nbrTri*9;
NumNonA=0;
Ai=zeros(NbrNon,1);
Aj=zeros(NbrNon,1);
Av=complex(zeros(NbrNon,1));
NumNonB=0;
Bi=zeros(NbrNon*4,1);
Bj=zeros(NbrNon*4,1);
Bv=complex(zeros(NbrNon*4,1));

PECNode=sort(unique([PECNode;CnzDstIndex]));
PECEdge=sort(unique([PECEdge;CntDstIndex]));
mesh.nbrPECIndex=length(PECNode);
PECNode=[PECNode;0];
PECEdge=[PECEdge;0];

%loop tri
for n=1:mesh.nbrTri
    %coordinate of tri's vertex
    x=zeros(3,1);y=zeros(3,1);l=ones(3,1);
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
    mu=physic.mur(domain);
    epsilon=physic.epsilonr(domain);

    %submatrix
    St=zeros(3,3);Tte=zeros(3,3);
    Sz=zeros(3,3);Tz=zeros(3,3);
    G=zeros(3,3);Ttu=zeros(3,3);
    for i=1:3
       for j=1:3
          for k=1:IntegralOrder
              St(i,j)=St(i,j)+pt(k)*DetJac/mu*dot(curlEt(i,:,k),curlEt(j,:,k));
              Tte(i,j)=Tte(i,j)+pt(k)*DetJac*k0*k0*epsilon*dot(Et(i,:,k),Et(j,:,k));
              Sz(i,j)=Sz(i,j)+pt(k)*DetJac/mu*dot(gradEz(i,:,k),gradEz(j,:,k));
              Tz(i,j)=Tz(i,j)+pt(k)*DetJac*k0*k0*epsilon*dot(Ez(i,:,k),Ez(j,:,k));
              G(i,j)=G(i,j)+pt(k)*DetJac/mu*dot(Et(i,:,k),gradEz(j,:,k));
              Ttu(i,j)=Ttu(i,j)+pt(k)*DetJac/mu*dot(Et(i,:,k),Et(j,:,k));
          end
       end
    end
    Gt=G';
    Ae=[zeros(3,3) zeros(3,3);zeros(3,3) St-Tte];
    Be=[Sz-Tz Gt;G Ttu];
    
    %treat rotation periodic boundary
    MappingIndexS=zeros(3,1);MappingIndexV=zeros(3,1);
    P=eye(6,6);Pinv=eye(6,6);
    for i=1:3
        MappingIndexS(i)=mesh.tri(n,i);
        MappingIndexV(i)=mesh.edgesOfTri(n,i);
    end
    for i=1:3
       if (isempty(find(CnzDstIndex==MappingIndexS(i), 1)))
           
       else
           temp=find(CnzDstIndex==MappingIndexS(i), 1);
           MappingIndexS(i)=CnzSrcIndex(temp);
           P(i,i)=exp(1i*physic.m*Lambda);
           Pinv(i,i)=exp(-1i*physic.m*Lambda);
       end
       if (isempty(find(CntDstIndex==MappingIndexV(i), 1)))
           
       else
           temp=find(CntDstIndex==MappingIndexV(i), 1);
           MappingIndexV(i)=CntSrcIndex(temp);
           P(i+3,i+3)=exp(1i*physic.m*Lambda)*CnSign(temp);
           Pinv(i+3,i+3)=exp(-1i*physic.m*Lambda)*CnSign(temp);
       end
    end
    Ae=Pinv*Ae*P;
    Be=Pinv*Be*P;
    
    %treat PEC
    for i=1:3
        for j=1:3
            MappingIndexSi=MappingIndexS(i);
            MappingIndexSj=MappingIndexS(j);
            MappingIndexVi=MappingIndexV(i);
            MappingIndexVj=MappingIndexV(j);
            numSi=1;numSj=1;numVi=1;numVj=1;
            flagSi=0;flagSj=0;flagVi=0;flagVj=0;

            for k=1:mesh.nbrPECIndex  %PECNode,PECEdge
                if MappingIndexSi>PECNode(k)
                    numSi=numSi+1;
                end
                if MappingIndexSj>PECNode(k)
                    numSj=numSj+1;
                end
                if MappingIndexVi>PECEdge(k)
                    numVi=numVi+1;
                end
                if MappingIndexVj>PECEdge(k)
                    numVj=numVj+1;
                end
            end

            if MappingIndexSi==PECNode(numSi)
                flagSi=1;
            end
            if MappingIndexSj==PECNode(numSj)
                flagSj=1;
            end
            if MappingIndexVi==PECEdge(numVi)
                flagVi=1;
            end
            if MappingIndexVj==PECEdge(numVj)
                flagVj=1;
            end

            MappingIndexSi=MappingIndexSi-numSi+1;
            MappingIndexSj=MappingIndexSj-numSj+1;
            MappingIndexVi=MappingIndexVi-numVi+1+mesh.nbrVertex-mesh.nbrPECIndex;
            MappingIndexVj=MappingIndexVj-numVj+1+mesh.nbrVertex-mesh.nbrPECIndex;
            

            if (flagVi+flagVj)==0
                NumNonA=NumNonA+1;
                Ai(NumNonA)=MappingIndexVi;
                Aj(NumNonA)=MappingIndexVj;
                Av(NumNonA)=Ae(i+3,j+3);
                
                NumNonB=NumNonB+1;
                Bi(NumNonB)=MappingIndexVi;
                Bj(NumNonB)=MappingIndexVj;
                Bv(NumNonB)=Be(i+3,j+3);
            end
            if (flagSi+flagSj)==0
                NumNonB=NumNonB+1;
                Bi(NumNonB)=MappingIndexSi;
                Bj(NumNonB)=MappingIndexSj;
                Bv(NumNonB)=Be(i,j);
            end
            if (flagVi+flagSj)==0
                NumNonB=NumNonB+1;
                Bi(NumNonB)=MappingIndexVi;
                Bj(NumNonB)=MappingIndexSj;
                Bv(NumNonB)=Be(i+3,j);
            end
            if (flagSi+flagVj)==0
                NumNonB=NumNonB+1;
                Bi(NumNonB)=MappingIndexSi;
                Bj(NumNonB)=MappingIndexVj;
                Bv(NumNonB)=Be(i,j+3);
            end
        end
    end
    
end

PECNode=PECNode(1:end-1,1);
PECEdge=PECEdge(1:end-1,1);

Ai=Ai(1:NumNonA);
Aj=Aj(1:NumNonA);
Av=Av(1:NumNonA);
Bi=Bi(1:NumNonB);
Bj=Bj(1:NumNonB);
Bv=Bv(1:NumNonB);
AA = sparse(Ai,Aj,Av);
BB = sparse(Bi,Bj,Bv);

A=BB;
B=BB+1/physic.theta*AA;

mesh.PECNode=PECNode;
mesh.PECEdge=PECEdge;
mesh.CnzDstIndex=CnzDstIndex;
mesh.CntDstIndex=CntDstIndex;
mesh.CnzSrcIndex=CnzSrcIndex;
mesh.CntSrcIndex=CntSrcIndex;
mesh.CnSign=CnSign;

solver.r=symrcm(B);
solver.B=B(solver.r,solver.r);
solver.A=A(solver.r,solver.r);

end