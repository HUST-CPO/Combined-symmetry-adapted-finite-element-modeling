clc
clear all
close all

addpath('function');
addpath('kernel');
addpath('mesh');
addpath('post');
addpath('solver');

%read mesh data
load HollowCoreFiber.mat
Mesh.vertex=Mesh.vertex*1e-6;
Mesh.nbrVertex=length(Mesh.vertex);
Mesh.nbrTri=length(Mesh.tri);
%treat mesh data
Mesh=GetEdge(Mesh);
%boundary condition setting
Mesh.PECIndex=FindIndex([17 18 63 64],Mesh);

%physic model
physic.c_const=299792458;
physic.lam0=1.55e-06;%wavelength
physic.k0=2*pi/physic.lam0;
%material
physic.PML=[10];
physic.PMLData=[5.95e-5  7.75e-6];%r0 dr
physic.epsilonr=ones(19,1);%Corresponding to triID
physic.epsilonr(2:11)=1.4378*1.4378;
physic.mur=ones(19,1);%Corresponding to triID
physic.theta=physic.k0*physic.k0*1;

%assemble
tic
[solver,Mesh]=Assemble_standardFEM(Mesh,physic);
solver=Assemble_standardFEM_PEC(solver,Mesh);
disp('Matrix Assemble:')
toc

%solution
tic;
solver.nbrMode=40;
solver=solver_standardFEM(solver,physic.theta,physic.k0);
disp('solution:')
toc;
solver.neff

%post
Ez=solver.x(1:Mesh.nbrVertex,:);
Et=solver.x(Mesh.nbrVertex+1:end,:);
Ex=zeros(Mesh.nbrVertex,solver.nbrMode);Ey=zeros(Mesh.nbrVertex,solver.nbrMode);
for i=1:solver.nbrMode
    [Ex(:,i),Ey(:,i)]=GetExEy(Mesh.vertex,Mesh.tri,Mesh.edgesOfTri,Et(:,i));
    Ez(:,i)=Ez(:,i)*solver.gamma(i)*1i;
end

%plot
nn=1;
PlotnormE(Mesh.vertex,Mesh.tri,Ex(:,nn),Ey(:,nn),Ez(:,nn),0);%figure2
% PlotE(Mesh.vertex,Mesh.tri,abs(Ez(:,nn)),0);%figure1 
% PlotTri(Mesh.vertex,Mesh.tri);%figure3