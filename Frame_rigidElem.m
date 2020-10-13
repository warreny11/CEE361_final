% Frame, 4 towers - CEE 361
% An attempt at creating four towers of the Fred Hartman
% Author(s) : Jessica Chen, Warren Yuan
% Date : 10/12 - 

clc;
clear;
close all;

% num of towers
twr = 2;    % for test

% Define global parameters for Frame (single tower)
nsd = 3;        % number of spatial dimensions 
nel = 6;        % number of elements
nnp = 6;        % number of nodes
ndf = 7;        % number of degrees of freedom (doesn't change, ever!)

% Nodal definitions
xn = zeros(nnp,nsd);	% xyz nodal coordinates [m]

% (right top) w/ rigid elems
xn = [381,12,-54.25; 
      381,24,0; 381,23,0; 381,1,0; 381,0,0;
      381,12,80];  
  
% Supported Dofs
idb = zeros(nnp,ndf);	% index of dofs - supported
idb(1,:) = 1;           % fixed support at node 1

% Element Connections
ien = zeros(nel,2); % index of element nodes
% for a single tower
ien = [1 2; 1 5; 2 3; 2 6; 4 5; 5 6];    % n-connection array

E = 50000;  % E modulus [N/mm^2]
A = 1e3;    % A area [mm^2]
I = 1e9;    % I moment [mm^4]

prop = zeros(nel,16);	% element properties
prop(:,1) = 3;	% element type [3 = frame]
prop(:,2) = E;	% E modulus [N/mm^2]

% the normal A and I 
prop(:,3) = A;    % A area [mm^2]
prop(:,4) = I;    % I moment [mm^4]

% rigid elem. A and I 
prop([3,4],3) = A*10;    % A area [mm^2]
prop([3,4],4) = I*10;    % I moment [mm^4]

if twr == 1 % for one tower, the base case
    disp("all good to go :)");
      
elseif twr == 2 % for two towers
    % right, two
    
    % nodes are double (what about the middle node? it's not a hinge!)
    % elements double
    nnp2 = 2*nnp;
    nel2 = 2*nel;
    
    % create y-translating vector
    y = zeros(nnp,nsd);
    y(:,2) = 24;
    
%     size(y)
%     size(xn)
    
    % translate and concatenate for second tower set of nodes
    xn2 = xn - y;
    xn = [xn; xn2];
    
    % init new ien for more towers
    ien = [nel2,2];
    ien = [1 2; 1 5; 2 3; 2 6; 4 5; 5 6; 5 8;  
           7 8; 7 11; 8 9; 8 12; 10 11; 11 12];
       
    prop = zeros(nel,16);	% element properties
    prop(:,1) = 3;	% element type [3 = frame]
    prop(:,2) = E;	% E modulus [N/mm^2]

    % the normal A and I 
    prop(:,3) = A;    % A area [mm^2]
    prop(:,4) = I;    % I moment [mm^4]
       
    % rigid elem. A and I 
    prop([5,8],3) = A*10;    % A area [mm^2]
    prop([5,8],4) = I*10;    % I moment [mm^4]
   
elseif twr == 4 % for all four towers
    % nodes are x 4
    % elements x 4
    nnp4 = 4*nnp;
    nel4 = 4*nel;
    
    % create y-translating vector
    y = zeros(nnp,nsd);
    y(:,2) = 24;
    
    % translate and concatenate for second tower set of nodes
    xn2 = xn - y;
    xn4 = [ -xn(:,1), xn(:,2:3);-xn2(:,1), xn2(:,2:3)]; 
    xn = [xn; xn2; x4];
    
    % init new ien for more towers
    ien = [nel4,2];
    ien = [1 2; 1 5; 2 3; 2 6; 4 5; 5 6; 5 8;  
           7 8; 7 11; 8 9; 8 12; 10 11; 11 12];
       
    prop = zeros(nel,16);	% element properties
    prop(:,1) = 3;	% element type [3 = frame]
    prop(:,2) = E;	% E modulus [N/mm^2]

    % the normal A and I 
    prop(:,3) = A;    % A area [mm^2]
    prop(:,4) = I;    % I moment [mm^4]
       
    % rigid elem. A and I 
    prop([5,8],3) = A*10;    % A area [mm^2]
    prop([5,8],4) = I*10;    % I moment [mm^4]
end



% writeDXF for this example
filenm = "frame_rigid_twr" + twr;
writeDXF(filenm,xn(:,1),xn(:,2),xn(:,3),ien);
