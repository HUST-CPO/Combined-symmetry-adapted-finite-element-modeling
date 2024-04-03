function NodeIndex=FindIndex(flag,mesh)

if isempty(flag)
    NodeIndex=[];
    return
end

index=find(mesh.BedgeID==flag(1));
if length(flag)>1
    for i=2:length(flag)
        index=[index;find(mesh.BedgeID==flag(i))];
    end
end

index=sort(index);
boundary=mesh.Bedge(index,:);
boundary=sort(boundary,2);

%vertex
NodeIndex=unique(sort([boundary(:,1);boundary(:,2)]));

%Bedge
[~,index]=ismember(boundary,mesh.edges,'rows');
index=index+mesh.nbrVertex;

NodeIndex=sort([NodeIndex ; index]);
