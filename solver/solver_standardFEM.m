function solver=solver_standardFEM(solver,theta,k0)

%Prioritize calculating the base mode
A=solver.B;
B=solver.B+1/theta*solver.A;

%Matrix reordering
r=symrcm(B);
A=A(r,r);
B=B(r,r);

%solution
tic
[x,GAMA_sq]=eigs(A,B,solver.nbrMode,'largestreal');
disp('solution:')
toc

%Recovery sort
q=zeros(length(r),1);
for i=1:length(r)
    q(r(i))=i;
end
for i=1:solver.nbrMode
   x(:,i)=x(q,i); 
end

%treat PEC
solver.x=solver.P*x;

%neff
lambda=diag(GAMA_sq);
beta=sqrt(theta-theta./lambda);
solver.neff=beta./k0;
solver.gamma=solver.neff.*k0;

end