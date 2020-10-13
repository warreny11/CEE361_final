% Frame, 4 towers - CEE 361
% An attempt at creating four towers of the Fred Hartman
% Author(s) : Jessica Chen, Warren Yuan
% Date : 10/12 - 

clc;
clear;
close all;

% num of towers
twr = 4;

% Define global parameters for Frame (single tower)
nsd = 3;        % number of spatial dimensions 
nel = 6;        % number of elements
nnp = 6;        % number of nodes
ndf = 7;        % number of degrees of freedom (doesn't change, ever!)

% init bridge dimensions
dhs = 381;      % distance to half span (dhs) [m]
rel = 1;        % rigid element length (rel) [m]
hod = -54.25;   % height of deck, where deck is at 0m (hod) [m] 
wod = 24;       % width of deck, maximum width of tower (wod) [m]
hot = 80;       % max height of tower (hot) [m]

% Nodal definitions
xn = zeros(nnp,nsd);	% xyz nodal coordinates [m]

% (right top) w/ rigid elems
xn = [dhs,wod/2,hod; 
      dhs,wod,0; dhs,wod-rel,0; dhs,rel,0; dhs,0,0;
      dhs,wod/2,hot];  
  
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
    nnp = 2*nnp;
    nel = 2*nel;
    
    % translate and concatenate for second tower set of nodes
    xn2 = xn;
    xn2(:,2) = -xn(:,2);
    xn = [xn; xn2];
    
    % init new ien for more towers
    ien = [nel,2];
    ien = [1 2; 1 5; 2 3; 2 6; 4 5; 5 6;  
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
    nnp = 4*nnp;
    nel = 4*nel;
    
    % translate and concatenate for second tower set of nodes
    % y translation
    xn2 = xn;
    xn2(:,2) = -xn(:,2);
    
    % concatenate translated tower
    xn = [xn; xn2];
    
    % x translation
    xn4 = xn; 
    xn4(:,1) = -xn(:,1);
    
    % concatenate translated tower
    xn = [xn; xn4];
    
    % init new ien for more towers
    ien1 = [1 2; 1 5; 2 3; 2 6; 4 5; 5 6];
    ien2 = ien1 + 6;
    ien3 = ien2 + 6;
    ien4 = ien3 + 6;
    
    ien = [nel,2];
    ien = [ien1;ien2;ien3;ien4];
       
    % this still needs fixing!
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
