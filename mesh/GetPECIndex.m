function [PECNode,PECEdge]=GetPECIndex(Boundary,Edges)

PECNode=[Boundary(:,1);Boundary(:,2)];
PECNode=unique(PECNode);
PECNode=sort(PECNode);


[~ ,PECEdge] = ismember(Boundary,Edges,'rows');
PECEdge=sort(PECEdge);

end