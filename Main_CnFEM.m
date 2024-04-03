clc
clear all
close all

addpath('function');
addpath('kernel');
addpath('mesh');
addpath('post');
addpath('solver');

%read mesh data
load OpticalFiber_model_CnFEM.mat
Mesh.nbrVertex=length(Mesh.vertex);
Mesh.nbrTri=length(Mesh.tri);
%treat mesh data
Mesh=GetEdge(Mesh);
%boundary condition setting
Mesh.PECID=[6];
Mesh.SrcID=[1 4];
Mesh.DstID=[2,5];
[Mesh.SrcBoundary,Mesh.DstBoundary,Mesh.OrigIndex]=GetPeriodicBoundIndex(Mesh);


%physic model
physic.m=0;% m=0,1,..N-1  category of mode
physic.N=6;%C6
physic.c_const=299792458;
physic.lam0=1.55E-6;%wavelength
physic.k0=2*pi/physic.lam0;
%material
physic.epsilonr=ones(2,1);%Corresponding to triID
physic.epsilonr(1)=2.1025;
physic.epsilonr(2)=1.96;
physic.mur=ones(2,1);%Corresponding to triID
physic.maxEps=max(physic.epsilonr);
physic.theta=physic.k0*physic.k0*physic.maxEps;

%assemble
tic
[solver,Mesh]=Assemble_CnFEM(Mesh,physic);
disp('Matrix Assemble:')
toc

%solution
solver.nbrMode=5;
solver=solver_CnFEM(solver,physic,Mesh);
solver.neff

%post
for i=1:solver.nbrMode
    Ez(:,i)=solver.x(1:Mesh.nbrVertex,i)*(1i*solver.gamma(i));
end
Et=solver.x(Mesh.nbrVertex+1:Mesh.Dof,:);

%plot
nn=1;
% plotEz(Nodes,Elements,abs(Ez(:,nn)),0);
[Ex,Ey]=GetExEy(Mesh.vertex,Mesh.tri,Mesh.edgesOfTri,Et(:,nn));
PlotE(Mesh.vertex,Mesh.tri,imag(Ez(:,nn)),0);%figure1 
PlotnormE(Mesh.vertex,Mesh.tri,Ex(:,1),Ey(:,1),Ez(:,nn),0);%figure2
PlotTri(Mesh.vertex,Mesh.tri);%figure3