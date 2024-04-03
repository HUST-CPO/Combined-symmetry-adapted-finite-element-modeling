function [Ex,Ey]=GetExEy(node,elem,EdgesOfElements,Et)

Elements = elem';
NbrElement=max(size(Elements));

u=[0,1,0];v=[0,0,1];
NbrNodes=length(node);
Ex=zeros(NbrNodes,1);
Ey=zeros(NbrNodes,1);
numEt=zeros(NbrNodes,1);
for n=1:NbrElement
    x=zeros(3,1);y=zeros(3,1);l=zeros(3,1);
    x(1)=node(Elements(1,n),1);y(1)=node(Elements(1,n),2);
    x(2)=node(Elements(2,n),1);y(2)=node(Elements(2,n),2);
    x(3)=node(Elements(3,n),1);y(3)=node(Elements(3,n),2);
    l(1)=sqrt((x(1)-x(2))*(x(1)-x(2))+(y(1)-y(2))*(y(1)-y(2)));
    l(2)=sqrt((x(1)-x(3))*(x(1)-x(3))+(y(1)-y(3))*(y(1)-y(3)));
    l(3)=sqrt((x(2)-x(3))*(x(2)-x(3))+(y(2)-y(3))*(y(2)-y(3)));
    Jac=zeros(3,3);
    Jac(1,1)=x(2)-x(1);Jac(1,2)=y(2)-y(1);
    Jac(2,1)=x(3)-x(1);Jac(2,2)=y(3)-y(1);Jac(3,3)=1;
    InvJac=inv(Jac);
    %BF
    temp=zeros(3,1);et=zeros(3,3,3);
    for i=1:3
       for j=1:3
           [temp(1),temp(2)]=BF_Et(j,u(i),v(i));temp(3)=0;
            et(i,:,j)=InvJac*temp*l(j);
       end
    end

    for i=1:3
        for j=1:3
            Ex(Elements(i,n))=Ex(Elements(i,n))+et(i,1,j)*Et(EdgesOfElements(n,j));
            Ey(Elements(i,n))=Ey(Elements(i,n))+et(i,2,j)*Et(EdgesOfElements(n,j));
        end
        numEt(Elements(i,n))=numEt(Elements(i,n))+1;
    end
end
Ex=Ex./numEt;
Ey=Ey./numEt;

end
