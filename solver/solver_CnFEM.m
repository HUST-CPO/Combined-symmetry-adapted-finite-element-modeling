function solver=solver_CnFEM(solver,physic,mesh)

%solution
tic
[x,GAMA_sq]=eigs(solver.A,solver.B,solver.nbrMode,'largestreal');
disp('solution:')
toc

%Recovery sort
q=zeros(length(solver.r),1);
for i=1:length(solver.r)
    q(solver.r(i))=i;
end
for i=1:solver.nbrMode
   x(:,i)=x(q,i); 
end

%neff
lambda=diag(GAMA_sq);
beta=sqrt(physic.theta-physic.theta./lambda);
solver.neff=beta./physic.k0;
solver.gamma=solver.neff.*physic.k0;

%treat boundary condition
for i=1:mesh.nbrPECIndex
    n=mesh.PECNode(i);
    x=[x(1:n-1,:);zeros(1,solver.nbrMode);x(n:end,:)];
end
for i=1:mesh.nbrPECIndex
    n=mesh.PECEdge(i)+mesh.nbrVertex;
    if n<=length(x)
        x=[x(1:n-1,:);zeros(1,solver.nbrMode);x(n:end,:)];
    else
        x=[x(1:n-1,:);zeros(1,solver.nbrMode);];
    end
end
for i=1:length(mesh.CnzDstIndex)
    n1=mesh.CnzDstIndex(i);
    n2=mesh.CnzSrcIndex(i);
    x(n1,:)=x(n2,:)*exp(1i*physic.m*2*pi/physic.N);  
end
for i=1:length(mesh.CntDstIndex)
    n1=mesh.CntDstIndex(i)+mesh.nbrVertex;
    n2=mesh.CntSrcIndex(i)+mesh.nbrVertex;
    x(n1,:)=x(n2,:)*mesh.CnSign(i)*exp(1i*physic.m*2*pi/physic.N);  
end

solver.x=x;
