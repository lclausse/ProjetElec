clear;
clc;

[txt_theta, txt_phi, txt_amplitude, txt_phase] = importfileDiag("Diagramme.txt",1,1000000);
% theta : 0 -> 180 par pas de 5 -> 37 elements
% phi : 0 -> 355 par pas de 5   -> 72 elements

data = reshape([txt_theta txt_phi txt_amplitude txt_phase],[37,72,4]);
% Array 2D -> max et min dans les deux directions
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


[amp, ph] = interpolate(13, 9)







% --------- Plot du mesh en cartesien ---------
% X : phi (horizontal), Y : theta (vertical)
%mesh(data(:,:,2), data(:,:,1), data(:,:,3));
% ---------------------------------------------


% --------- Plot du mesh en spherique ---------
% X : phi, Y : theta
[x,y,z] = sph2cart(deg2rad(data(:,:,2)), deg2rad(data(:,:,1)), data(:,:,3));

% Uncomment for amplitude
%{
subplot(1,2,1);
C = data(:,:,4);
caxis([240 300])
surf(x,y,z,C,'FaceAlpha',1)
colorbar

subplot(1,2,2);
plot3(x,y,z)
%}


% Uncomment for phase
%{
subplot(1,2,1);
mesh(data(:,:,2), data(:,:,1), data(:,:,4));

subplot(1,2,2);
plot(data(:,:,4))
%}
% ---------------------------------------------

%Positions d'un récepteur et d'un transmetteur:
% [xt, yt, zt]
%{
pos_t = [1, 1, 1];
pos_r = [2, 2.5, 3];

scatter3(pos_t(1), pos_t(2), pos_t(3), 'filled', 'r');
hold on;
scatter3(pos_r(1), pos_r(2), pos_r(3), 'filled');
hold on;
axis([0.5 2.5 0.5 3.5 0 3.5]);
%}

%{
plot3([pos_t(1) pos_r(1)], [pos_t(2) pos_r(2)], [pos_t(3) pos_r(3)]);
hold on;

% Repère
l = 1;
plot3([pos_t(1) pos_t(1)+l], [pos_t(2) pos_t(2)], [pos_t(3) pos_t(3)],'r');
text(pos_t(1)+ 1.1*l, pos_t(2), pos_t(3), 'x')
hold on;
plot3([pos_t(1) pos_t(1)], [pos_t(2) pos_t(2)+l], [pos_t(3) pos_t(3)],'r');
text(pos_t(1), pos_t(2)+ 1.1*l, pos_t(3), 'y')
hold on;
plot3([pos_t(1) pos_t(1)], [pos_t(2) pos_t(2)], [pos_t(3) pos_t(3)+l],'r');
text(pos_t(1), pos_t(2), pos_t(3)+ 1.1*l, 'z')
hold on;

plot3([pos_t(1) pos_r(1)], [pos_t(2) pos_r(2)], [pos_r(3) pos_r(3)],'k--');
hold on;
plot3([pos_t(1) pos_t(1)], [pos_t(2) pos_t(2)], [pos_t(3) pos_r(3)],'k--');
hold on;
plot3([pos_t(1) pos_r(1)], [pos_t(2) pos_r(2)], [pos_t(3) pos_t(3)],'k--');
hold on;
plot3([pos_r(1) pos_r(1)], [pos_r(2) pos_r(2)], [pos_t(3) pos_r(3)],'k--');
hold on;
plot3([pos_t(1) pos_r(1)], [pos_t(2) pos_t(2)], [pos_t(3) pos_t(3)],'k--');
hold on;
plot3([pos_r(1) pos_r(1)], [pos_t(2) pos_r(2)], [pos_t(3) pos_t(3)],'k--');

retrait = 0.2;
text(pos_t(1) - retrait, pos_t(2) - retrait, pos_t(3) - 2*retrait, '(x_t, y_t, z_t)')
text(pos_r(1) + 0.05, pos_r(2) + 0.05, pos_r(3) + retrait, '(x_r, y_r, z_r)')
%}

% Retourne l'amplitude et la phase après interpolation
function [amp, pha] = interpolate(theRand, phiRand)
    index_theRand = (theRand / 5) + 1;
    index_phiRand = (phiRand / 5) + 1;
    
    index_thei1 = floor(index_theRand);
    index_thei2 = ceil(index_theRand);
    index_phii1 = floor(index_phiRand);
    index_phii2 = ceil(index_phiRand);
    
    
    
    amp = 0;
    pha = 0;
end

