function mesh=GetEdge(mesh)

el2no=mesh.tri';
n1=el2no([1 1 2],:);
n2=el2no([2 3 3],:);
el_ed2no_array=[n1(:) n2(:)];
[mesh.edges,~,mesh.edgesOfTri]=unique(el_ed2no_array,'rows');
mesh.nbrEdges=length(mesh.edges);
mesh.edgesOfTri=reshape(mesh.edgesOfTri,3,mesh.nbrTri)';

end