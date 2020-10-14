% Generate Mesh Deck - CEE 361
% A desperate attempt to simplify the mesh of the deck (for example,
% cylindrical shape is not needed.
% Author(s): Jessica Chen, Warren Yuan
% Date: 10/14 - 

% function genDeck takes width (W), length (L), num of X-nodes (X)
% and num of Y-nodes (Y). It outputs xn, ien.
function [xn,ien] = genMeshDeck(W,L,neX,neY)

nsd = 3;                % number of spatial dimensions (3D, so nsd = 3)
nel = neX*neY;          % number of elements 
nen = 4;  % number of element nodes (doesn't change, for square elements)
nnp = (neX+1)*(neY+1);	% number of nodal points

% nodal definitions
Xinc = W/neX;   % the x-increments (width/num of x elems)
Yinc = L/neY;   % the y-increments (length/num of y elems)

xn = zeros(nnp,nsd)   % position of nodes in 3D

% populate xn with nodal coordinates in 3D space
for i = 1: neX+1
  for j = 1: neY+1
    n = i + (j-1)*(neX+1);  % generates index for xn
    thetaN = Xinc*(i-1);    % generates 
    xn(n,:) = [(j-1)*Yinc-190.5 thetaN 12];
  end
end

% element definitions
ien = zeros(nel,nen);	% index of element nodes
for i = 1:neX
  for j = 1:neY
    e = i+(j-1)*neX;
    n1 = i+(j-1)*(neX+1); n2 = n1 + 1;
    n3 = i+j*(neX+1); n4 = n3 + 1;
    ien(e,:)     = [n1 n2 n4 n3];
  end
end

%xnew = isoProject([xn(:,1) xn(:,2) xn(:,3)]);
%writeDXF(strcat('deformed\test',num2str(neR)),xn(:,1),xn(:,2),xn(:,3),ien)
%writeDXF(strcat('deformed\test',num2str(neR)),xnew(:,1),xnew(:,2),xnew(:,1)*0,ien)
