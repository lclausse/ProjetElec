clear;
clc;

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

A1 = [1,0,5];
A2 = [1,3,5];