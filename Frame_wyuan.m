% Frame Geometry - CEE 361
% This file models a single frame for one of the 
% Fred Hartman Bridge towers.
% Authors : Jessica Chen and Warren Yuan
% Date : 10/9 - 

clc;
close all; 
clear;

% num of towers
twr = 2;

% Define global parameters 
nsd = 3; % number of spatial dimensions 
nel = 5*2; % number of elements
nnp = 4*2-1; % number of nodes
ndf = 7; % number of degrees of freedom (doesn't change, ever!)

% Nodal definitions
xn = zeros(nnp,nsd);	% xyz nodal coordinates [m]

if twr == 1
    % for one tower
    xn = [0,0,0;0,13,64.5;0,-13,64.5;0,0,129];    % nodes in 3D [m]
else
    % for two towers
    xn = [0,-13,0;0,13,0;0,-26,64.5;0,0,64.5;0,26,64.5;0,-13,129;0,13,129];
end

idb = zeros(nnp,ndf);	% index of dofs - supported
idb(1,[1,2]) = 1;  
idb(4,[1,6]) = 1;
% fixed support at nodes 1,4; pin @ n1, moment roller @ n4

ds = zeros(nnp,ndf);	% prescribed displacements at supports [m]
% no prescribed displacements

Pu = zeros(nnp,ndf); 	% applied forces (P) at unrestrained dofs [N]
P = -9*1000;    % concentrated load of 9 [kN]
Pu([2,3],2) = P; % nodal applied loads @ n2,3

% Element definitions
ien = zeros(nel,2); % index of element nodes
% for a single tower
% ien = [1 2; 1 3; 2 3; 2 4; 3 4];    % n-connection array
% for two towers
ien = [1 3; 1 4; 2 4; 2 5; 3 4; 3 6; 4 5; 4 6; 4 7; 5 7]; 

E = 50000;  % E modulus [N/mm^2]
A = 1e3;    % A area [mm^2]
I = 1e9;    % I moment [mm^4]

% A = A/10;   % A area, divided by 10
% I = I/10;   % I moment, divided by 10

prop = zeros(nel,16);	% element properties
prop(:,1) = 3;	% element type [3 = frame]
prop(:,2) = E;	% E modulus [N/mm^2]
prop(:,3) = A;    % A area [mm^2]
prop(:,4) = I;    % I moment [mm^4]
  
% writeDXF for this example
filenm = twr + "_" + nnp + "_" + nel + "_";
writeDXF(filenm,xn(:,1),xn(:,2),xn(:,3),ien);

% % RUN ANALYSIS
% [results,process] = runAnalysis(Pu,ds,xn,prop,idb,ien);
% [F,Rs,Fe,Fi,d,du,de] = deal(results{:});
% [Kuu,Ke,ke,Te,ied,idu,ids] = deal(process{:});