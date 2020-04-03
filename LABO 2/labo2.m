clear;
clc;

load("data_labo_reflexion.mat");

%subplot(2, 2, 1);
plot(Reference);
%subplot(2, 2, 2);
%plot(Corr7);



function [vec] = fourier(r)
    L = length(r);
    R = fftshift(fft(r));
    vec = abs(R/L);
end