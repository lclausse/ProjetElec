clear;
clc;

%Q3.3.2 a)

t = 0:10e-3:1;
Fs=10e8;

% Fréquences acceptées pour le sinc
f = linspace(-3*10e9, 3*10e9,10000);
fOk = [-0.8*10e9 0.8*10e9];
% Amplitude de la pulsation
R = zeros(1,length(f));
R(find(f<fOk(2) & f>fOk(1))) = 10^-5;

subplot(2,1,1);
plot(f,R);

r = real(fftshift(ifft(R)));
subplot(2,1,2);
plot(r);

%Q3.3.2 b)

A1 = [1,0,5];
A2 = [1,3,5];
PI = [1,3/2,0];



theta = (180*atan(3/10))/pi;

eperp = cross(PI-A1,PI-A2)/norm(cross(PI-A1,PI-A2));
eipara = cross(eperp,A1-PI)/norm(cross(eperp,A1-PI));
erpara = cross(eperp,PI-A2)/norm(cross(eperp,PI-A2));

Rpara = (3.3*cos(theta)-sqrt(3.3-sin(theta)^2))/(3.3*cos(theta)+sqrt(3.3-sin(theta)^2));
Rperp = (cos(theta)-sqrt(3.3-sin(theta)^2))/(cos(theta)+sqrt(3.3-sin(theta)^2));
Einc = [0,0,1];
Erec = Einc;

E=Rpara*dot(Einc,eipara)*dot(erpara,Erec)+Rperp*dot(Einc,eperp)*dot(eperp,Erec)
