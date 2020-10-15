% Complete Bridge - CEE 361 
% This file combines the dxf/geometric outputs of several files (endspans
% midspan, and the frame, to obtain a single, combined dxf.
% Author(s): Jessica Chen, Warren Yuan
% Date: 10/14 - 

clc;
clear;
close all;

% file marker num
i = 11;

% for combining decks
neR_end = 4;
neY_end = 12;
L_end = 180;
W_end = 22;

neR_mid = 4;
neY_mid = 24;
L_mid = 381;
W_mid = 22;

[xn1,ien1] = Deck_endspan(neR_end,neY_end,L_end,W_end);   
% first endspan

[xn2,ien2] = Deck_midspan(neR_mid,neY_mid,L_mid,W_mid);    
% midspan

xn3 = xn1;

xn1(:,1) = xn1(:,1) - L_end;   % shift over L_end
xn3(:,1) = -xn1(:,1);   % perform that y-axis inversion
ien3 = ien1;    % last endspan

[xn4, ien4] = Frame_tower; % towers (all four)

xn = [xn1+10; xn2+10; xn3+10];

offset1 = max(max(ien1));   % takes the offsets to ensure no overlap
offset2 = offset1 + max(max(ien2)); % adds subsequent offset
offset3 = offset2 + max(max(ien3)); % ""

ien2 = ien2 + offset1;
ien3 = ien3 + offset2;
ien4 = ien4 + offset3;  

ien4_zeros = zeros(length(ien4),2); %   still not working-need to figure out how to couple with deck.
ien4 = [ien4, ien4_zeros];
    
ien = [ien1; ien2; ien3];
% ien = [ien1; ien2+offset1; ien3+offset2];

filenm = "completeBridge_" + i + "th_run";
writeDXF(filenm,xn(:,1),xn(:,2),xn(:,3),ien);
