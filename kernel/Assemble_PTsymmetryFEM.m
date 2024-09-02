function solver=Assemble_PTsymmetryFEM(solver,mesh)

nn=length(solver.A);
Ar=real(solver.A);
Ai=imag(solver.A);
Br=real(solver.B);
Bi=imag(solver.B);
A=sparse([Ar -Ai;Ai Ar]);
B=sparse([Br -Bi;Bi Br]);

% delelteIndex=unique(sort([mesh.PECIndex;mesh.PECIndex+nn;mesh.PTIndex+nn]));
delelteIndex=unique(sort([mesh.PECIndex;mesh.PECIndex+nn;mesh.PTIndex]));
Index=1:(nn*2);
Index(delelteIndex)=[];

P=speye(nn*2);
P=P(:,Index);

solver.P=P;
solver.A=P'*A*P;
solver.B=P'*B*P;

end