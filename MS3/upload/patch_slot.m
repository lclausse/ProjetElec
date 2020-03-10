function [PHI_new,THT_new,R_new] = patch_slot(PHI,THT,R)
% phi must contain both 0 and pi such that the plots do not have a slot
PHI_new=[PHI PHI(:,1)+pi];
THT_new=[THT THT(:,1)];
R_new=[R flipud(R(:,1))];

