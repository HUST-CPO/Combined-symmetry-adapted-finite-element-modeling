clc
clear all

%loading COMSOL project
fileName="HollowCoreFiber";
model=mphload(fileName+".mph");

%set mesh size and perform mesh division
[~, m2]=mphmeshstats(model);

%process mesh data
Mesh.vertex=m2.vertex';%coordinate of nodes  unit:m
tri=m2.elem{2}+1;
Mesh.tri = sort(tri)';%tri
Mesh.triID=m2.elementity{2};%face ID of tri
Mesh.Bedge=m2.elem{1}'+1;%boundary edge
Mesh.BedgeID=m2.elementity{1};%edge Id of boundary edge

save(fileName,'Mesh');
