% Deck Geometry - CEE 361
% This code simulates a deck for the Fred-Hartman with four lanes and 
% twelve nodal connections. 
% Authors : Jessica Chen and Warren Yuan
% Date : 10/9 - 

clc;
close all; 
clear;

% 0. Convergence inputs
neR = 4;   % number of elements along x-axis (for each lane)
neY = 12;  % number of elements along y-axis (for each cable connection)
npl = 4; npm = 4; % membrane (lambda, mu) integration pts WHAT R THOSE?
npb = 4; npv = 1; % plate (bending, shear) integration pts

% 1. Global definitions :: these need to remain constant!
nsd = 3;	% number of spatial dimensions
ndf = 7;	% number of degrees of freedom

% Given constants?need to find our own
L = 381;    % length of deck over span (381 [m])
R = 10*L;   % radius of curvature, in this case, set arbitrarily large
t = 0.25;   % need to make sure that this/something is 28m total?
E = 4.32e8; 
v = 0; 
theta = 0.1*(pi/180);   % see textbook (p267), angle of arc 
q = -90; 
Pe = q*(L/neY/2)*(2*R*sin(1*pi/9/neR));

% 2. Nodal & 3. Element definitions
[xn,ien] = genMeshRoof(R,L/2,theta,neR,neY);
nnp = size(xn,1); % number of nodal points
nel = size(ien,1);  % number of elements
idb = zeros(nnp,ndf);	% index of dofs - supported
ds = zeros(nnp,ndf);	% prescribed displacements at supports (mm)
Pu = zeros(nnp,ndf); 	% applied forces (P) at unrestrained dofs (N)
prop = repmat([7 E 0 0 0 0 v t 0 0 0 0 npl npm npb npv],[nel 1]); % property matrix

for n = 1:nnp   % Set essential BCs
  if xn(n,2) == 0,   idb(n,[1 3 5]) = 1; end % simple-support
  if xn(n,1) == 0,   idb(n,[1 5 6]) = 1; end
  if xn(n,2) == L/2, idb(n,[2 4 6]) = 1; end
end
for e = 1:nel   % Set  natural BCs
  for i = 1:4
    Pu(ien(e,i),3) = Pu(ien(e,i),3) + Pe/4;
  end
end

% 4. writeDXF for this example
filenm = neR + "_" + neY + "_" + R + "_" + theta;
writeDXF(filenm,xn(:,1),xn(:,2),xn(:,3),ien);

% % 5. RUN ANALYSIS
% [results,process] = runAnalysis(Pu,ds,xn,prop,idb,ien);
% [F,Rs,Fe,Fi,d,du,de] = deal(results{:});
% [Kuu,Ke,ke,Te,ied,idu,ids] = deal(process{:});
% 
% % 6. Additional Post-Processing

