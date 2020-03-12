clear;
clc;

[txt_theta, txt_phi, txt_amplitude, txt_phase] = importfileDiag("Diagramme.txt",1,1000000);
% theta : 0 -> 180 par pas de 5 -> 37 elements
% phi : 0 -> 355 par pas de 5   -> 72 elements

data = reshape([txt_theta txt_phi txt_amplitude txt_phase],[37,72,4]);
% Array 2D -> max et min dans les deux directions
%plot(data(:,:,4));
phase_max = max(max(data(:,:,4)));
phase_min = min(min(data(:,:,4)));


% ------------ Raccord du diagramma -----------
% Pour faire le raccord, on ajoute une colonne en 360 degres
% phi : 0 -> 360 par pas de 5   -> 73 elements
data = cat(2, data, data(:,1,:));
% Il faut changer la valeur de phi au bout. Comme on a clone celle du
% debut, phi = 0 pour le moment -> pas bon
data(:,end,2) = 360;
% ---------------------------------------------


% --------- Plot du mesh en cartesien ---------
% X : phi (horizontal), Y : theta (vertical)
%mesh(data(:,:,2), data(:,:,1), data(:,:,3));
% ---------------------------------------------


% --------- Plot du mesh en spherique ---------
% X : phi, Y : theta
[x,y,z] = sph2cart(deg2rad(data(:,:,2)), deg2rad(data(:,:,1)), data(:,:,3));
%mesh(x,y,z)
%plot3(x,y,z)
C = data(:,:,4);
surf(x,y,z,C,'FaceAlpha',0.7)
colorbar
caxis([240 300])
% ---------------------------------------------

