function [SrcBoundary,DstBoundary,OrigIndex]=GetPeriodicBoundIndex(mesh)

SrcBoundaryFlag=mesh.SrcID;
DstBoundaryFlag=mesh.DstID;
%src
index=find(mesh.BedgeID==SrcBoundaryFlag(1));
for i=2:length(SrcBoundaryFlag)
    index=[index;find(mesh.BedgeID==SrcBoundaryFlag(i))];
end
index=sort(index);
SrcBoundary=mesh.Bedge(index,:);
SrcBoundary=sort(SrcBoundary,2);

%dst
index=find(mesh.BedgeID==DstBoundaryFlag(1));
for i=2:length(DstBoundaryFlag)
    index=[index;find(mesh.BedgeID==DstBoundaryFlag(i))];
end
index=sort(index);
DstBoundary=mesh.Bedge(index,:);
DstBoundary=sort(DstBoundary,2);

%origin
temp1=mesh.vertex(:,1)+mesh.vertex(:,2);
temp2=mesh.vertex(:,1).*mesh.vertex(:,2);
OrigIndex=find((temp1==0)&(temp2==0));