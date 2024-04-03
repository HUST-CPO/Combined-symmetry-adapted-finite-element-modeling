clc
clear all
close all

addpath('function');
addpath('kernel');
addpath('mesh');
addpath('post');
addpath('solver');

%read mesh data
load OpticalFiber_model_CnvFEM.mat
Mesh.nbrVertex=length(Mesh.vertex);
Mesh.nbrTri=length(Mesh.tri);
%treat mesh data
Mesh=GetEdge(Mesh); 
%boundary condition setting
% Mesh.PECIndex=findIndex([1,2,4,5,6],Mesh);%PEC-PEC
% Mesh.PECIndex=findIndex([6],Mesh);%PMC-PMC
% Mesh.PECIndex=findIndex([1,4,6],Mesh);%PEC-PMC
Mesh.PECIndex=FindIndex([2,5,6],Mesh);%PMC-PEC

%physic model
physic.c_const=299792458;
physic.lam0=1.55e-6;%wavelength
physic.k0=2*pi/physic.lam0;
physic.N=6;%C6v point group
physic.phi=2*pi/physic.N;
physic.m=0;%m=0 m=N/2
%material
physic.epsilonr=ones(2,1);%Corresponding to triID
physic.epsilonr(1)=2.1025;
physic.epsilonr(2)=1.96;
physic.mur=ones(2,1);%Corresponding to triID
physic.maxEps=max(physic.epsilonr);
physic.theta=physic.k0*physic.k0*physic.maxEps;

%assemble
tic
[solver,Mesh]=Assemble_standardFEM(Mesh,physic);
solver=Assemble_standardFEM_PEC(solver,Mesh);
disp('Matrix Assemble:')
toc

%solution
solver.nbrMode=3;
solver=solver_standardFEM(solver,physic.theta,physic.k0);
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
PlotE(Mesh.vertex,Mesh.tri,imag(Ez(:,nn)),0);%figure1 
PlotnormE(Mesh.vertex,Mesh.tri,Ex(:,nn),Ey(:,nn),Ez(:,nn),0);%figure2
PlotTri(Mesh.vertex,Mesh.tri);%figure3