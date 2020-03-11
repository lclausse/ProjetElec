clear;
clc;

[theta, phi, amplitude, phase] = importfileDiag("Diagramme.txt",1,1000000);
% theta : 0 -> 180 par pas de 5
% phi : 0 -> 360 par pas de 5

[phi, theta] = meshgrid(phi, theta);

%{
X = amplitude .* sin(phi) .* cos(theta);
Y = amplitude .* sin(phi) .* sin(theta);
Z = amplitude .* cos(phi);
%}

[x,y,z] = sph2cart(deg2rad(phi),deg2rad(theta),amplitude);



%scatter3(X,Y,Z,'.');
%warp(X,Y,Z,amplitude);
%surf(phi,theta,zeros(size(phi)));
surf(x,y,z);
%scatter(x,z,'.');
%plot3(x,y,z);
%Vq = interp3(X,Y,Z,V,Xq,Yq,Zq);

