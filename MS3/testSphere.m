clear;
clc;

[theta, phi, amplitude, phase] = importfileDiag("Diagramme.txt",1,1000000);

%scatter();


%[theta, phi] = meshgrid(theta, phi);

X = amplitude .* sin(phi) .* cos(theta(1));
Y = amplitude .* sin(phi) .* sin(theta(1));
Z = amplitude .* cos(phi);



scatter3(X,Y,Z,'.');


%surf(X,Y,Z);

%plot3(X,Y,Z);

%Vq = interp3(X,Y,Z,V,Xq,Yq,Zq);

%scatter3(X,Y,Z,'.');
