clc
clear all
close all

addpath('function');
addpath('kernel');
addpath('mesh');
addpath('post');
addpath('solver');

%read mesh data
load WaveguideWithPT-unit.mat
% load MT-symmtry.mat
Mesh.nbrVertex=length(Mesh.vertex);
Mesh.nbrTri=length(Mesh.tri);
%treat mesh data
Mesh=GetEdge(Mesh);
%boundary condition setting
Mesh.PTIndex=FindIndex([1],Mesh);
Mesh.PECIndex=FindIndex([3 9],Mesh);

%physic model
physic.c_const=299792458;
physic.lam0=1;%wavelength
physic.k0=2*pi/physic.lam0;
%material
physic.PML=[];
physic.epsilonr=ones(2,1);%Corresponding to triID
physic.epsilonr(1)=1;
physic.epsilonr(2)=10+0.4i;
physic.mur=ones(2,1);%Corresponding to triID
physic.theta=physic.k0*physic.k0*2.4*2.4 ;

%assemble
tic
[solver,Mesh]=Assemble_standardFEM(Mesh,physic);
solver=Assemble_PTsymmetryFEM(solver,Mesh);
disp('Matrix Assemble:')
toc

%solution
solver.nbrMode=2;
solver=solver_standardFEM(solver,physic.theta,physic.k0);
solver.neff

%post1
Ez=solver.x(Mesh.Dof+1:Mesh.Dof+Mesh.nbrVertex,:)-solver.x(1:Mesh.nbrVertex,:);
Et=solver.x(Mesh.Dof+Mesh.nbrVertex+1:Mesh.Dof+Mesh.Dof,:)-solver.x(Mesh.nbrVertex+1:Mesh.Dof,:);
Ex=zeros(Mesh.nbrVertex,solver.nbrMode);Ey=zeros(Mesh.nbrVertex,solver.nbrMode);
for i=1:solver.nbrMode
    [Ex(:,i),Ey(:,i)]=GetExEy(Mesh.vertex,Mesh.tri,Mesh.edgesOfTri,Et(:,i));
    Ez(:,i)=Ez(:,i)*solver.gamma(i)*1i;
end
%plot1
nn=1;
PlotnormE2(Mesh.vertex,Mesh.tri,Ex(:,nn),Ey(:,nn),Ez(:,nn),0);%figure2


%post2
Ez=solver.x(Mesh.Dof+1:Mesh.Dof+Mesh.nbrVertex,:)+solver.x(1:Mesh.nbrVertex,:);
Et=solver.x(Mesh.Dof+Mesh.nbrVertex+1:Mesh.Dof+Mesh.Dof,:)+solver.x(Mesh.nbrVertex+1:Mesh.Dof,:);
Ex=zeros(Mesh.nbrVertex,solver.nbrMode);Ey=zeros(Mesh.nbrVertex,solver.nbrMode);
for i=1:solver.nbrMode
    [Ex(:,i),Ey(:,i)]=GetExEy(Mesh.vertex,Mesh.tri,Mesh.edgesOfTri,Et(:,i));
    Ez(:,i)=Ez(:,i)*solver.gamma(i)*1i;
end
%plot2
nn=2;
PlotnormE2(-Mesh.vertex,Mesh.tri,Ex(:,nn),Ey(:,nn),Ez(:,nn),0);%figure2