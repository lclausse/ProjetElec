clear;
clc;

[txt_theta, txt_phi, txt_amplitude, txt_phase] = importfileDiag("Diagramme.txt",1,1000000);
% theta : 0 -> 180 par pas de 5 -> 37 elements
% phi : 0 -> 355 par pas de 5   -> 72 elements

global data;
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


[amp, ph] = interpolate(13, 9);







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
    global data;
    index_theRand = (theRand / 5) + 1;
    index_phiRand = (phiRand / 5) + 1;
    
    index_thei1 = floor(index_theRand);
    index_thei2 = ceil(index_theRand);
    index_phii1 = floor(index_phiRand);
    index_phii2 = ceil(index_phiRand);
    
    data1 = data(index_thei1, index_phii1,:);
    data2 = data(index_thei2, index_phii1,:);
    data3 = data(index_thei2, index_phii2,:);
    data4 = data(index_thei1, index_phii2,:);
    
    t1 = data1(1); t2 = data3(1);
    p1 = data1(2); p2 = data3(2);
    
    % Système linéaire de type A*x = B
    % -----------------------------------
    %{
    A = [t1 p1 t1*p1 1;
         t1 p2 t1*p2 1;
         t2 p1 t2*p1 1;
         t2 p2 t2*p2 1];
    B = [data1(3);
         data2(3);
         data3(3);
         data4(3)];
    coeff = A\B;    % Comme la matrice inverse, mais en plus efficace
    [thetaPlot,phiPlot] = meshgrid(t1:0.05:t2, p1:0.05:p2);
    ampPlot = coeff(1).*thetaPlot + coeff(2).*phiPlot + coeff(3).*thetaPlot.*phiPlot + coeff(4);
    %mesh(thetaPlot,phiPlot,ampPlot)
    %hold on;
    %}
    
    
    % http://supercomputingblog.com/graphics/coding-bilinear-interpolation/
    thetaAverage1Amp = ((p2-phiRand)/(p2-p1))*data1(3) + ((phiRand-p1)/(p2-p1))*data4(3);
    thetaAverage2Amp = ((p2-phiRand)/(p2-p1))*data2(3) + ((phiRand-p1)/(p2-p1))*data3(3);
    amp = ((t2-theRand)/(t2-t1))*thetaAverage1Amp + ((theRand-t1)/(t2-t1))*thetaAverage2Amp;
    
    thetaAverage1Phase = ((p2-phiRand)/(p2-p1))*data1(4) + ((phiRand-p1)/(p2-p1))*data4(4);
    thetaAverage2Phase = ((p2-phiRand)/(p2-p1))*data2(4) + ((phiRand-p1)/(p2-p1))*data3(4);
    pha = ((t2-theRand)/(t2-t1))*thetaAverage1Phase + ((theRand-t1)/(t2-t1))*thetaAverage2Phase;
    

    % ------------------- Plot -------------------
    scatter3(data1(1), data1(2), data1(3), 'filled', 'r');
    hold on;
    scatter3(data2(1), data2(2), data2(3), 'filled', 'r');
    hold on;
    scatter3(data3(1), data3(2), data3(3), 'filled', 'r');
    hold on;
    scatter3(data4(1), data4(2), data4(3), 'filled', 'r');
    hold on;
    scatter3(t1, phiRand, thetaAverage1Amp);
    hold on;
    scatter3(t2, phiRand, thetaAverage2Amp);
    hold on;
    scatter3(theRand, phiRand, amp);
    hold on;
    axis([data1(1)-1 data3(1)+1 data1(2)-1 data3(2)+1 0 20]);

end

