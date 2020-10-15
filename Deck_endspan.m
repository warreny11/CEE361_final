% Deck, Endspans - CEE 361
% The code for the endspans (single lane) of deck
% Author(s) : Jessica Chen, Warren Yuan
% Date : 10/12 - 

function [xn,ien] = Deck_endspan(neR,neY,L,W);
    % number of decks
    nod = 2;

    % 0. Convergence inputs
    npl = 2; npm = 2; % membrane (lambda, mu) integration pts
    npb = 1; npv = 1; % plate (bending, shear) integration pts

    % 1. Global definitions
    nsd = 3;	% number of spatial dimensions
    ndf = 7;	% number of degrees of freedom
    t = 0.304; 
    E = 4.32e8; 
    v = 0; % given values in structure
    q = -90; 
    Pe = q*(L/neY/2)*(2*W*sin(1*pi/9/neR));

    % 2. Nodal & 3. Element definitions
    if nod == 1
        [xn,ien] = genMeshDeck(W,L,neR,neY);
    end

    if nod == 2
        [xn,ien] = genMeshDeck(W,L,neR,neY);
        [xn2,ien2] = genMeshDeck(-W,L,neR,neY);
        xn(:,2) = [xn(:,2)+1];  % shifts over right deck to right 1m
        xn2(:,2) = [xn2(:,2)-1]; % shifts over left deck to left 1m
        xn = [xn;xn2];
        ien = [ien;ien2+max(max(ien))];
    end

    nnp = size(xn,1);   % number of nodal points
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

    xn(:,3) = xn(:,3) - W;

    filenm = "deck_endspan_" + nod + "dx1";
    writeDXF(filenm,xn(:,1),xn(:,2),xn(:,3),ien);
