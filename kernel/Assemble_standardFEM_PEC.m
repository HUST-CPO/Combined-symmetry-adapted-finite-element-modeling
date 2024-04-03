function solver=Assemble_standardFEM_PEC(solver,mesh)

solver.P=speye(mesh.Dof);

NbrPEC=length(mesh.PECIndex);
if NbrPEC>0
    for i=1:NbrPEC
        tempN=NbrPEC-i+1;
        tempN=mesh.PECIndex(tempN);
        solver.P(:,tempN)=[];
    end
end

solver.A=solver.P'*solver.A*solver.P;
solver.B=solver.P'*solver.B*solver.P;

end