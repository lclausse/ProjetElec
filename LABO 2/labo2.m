clear;
clc;

load("data_labo_reflexion.mat");

%subplot(2, 2, 1);
plot(Reference);
%subplot(2, 2, 2);
%plot(Corr7);






function [vec] = fourier(r)
    L = length(r);
    R = fft(r);
    P2 = abs(R/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    vec = P1;
end