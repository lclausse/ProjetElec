clear;
clc;

load("data_labo_reflexion.mat");

subplot(2, 2, 1);
plot(abs(Corr7));
decalsupp = 85;
attenuation = 0.66;

% reflet1 = bulbizarre(Reference,115,attenuation);
% reflet2 = bulbizarre(Reference,115,attenuation);
% reflet3 = bulbizarre(Reference,115,attenuation);
test1 = racaillou(113,attenuation,Reference);
subplot(2,2,2);
plot(abs(test1));
test2 = racaillou(115,attenuation,Reference);
subplot(2,2,3);
plot(abs(test2));
test3 = racaillou(118,attenuation,Reference); %118 ca donne bien
subplot(2,2,4);
plot(abs(test3));






function [vec] = racaillou(decalreflet,attenuation,Ref)
    L = length(Ref);
    removeEnd = L - decalreflet+1;
    RefRaccourci = Ref;
    RefRaccourci(removeEnd:end) = [];
    decalageReflet = zeros(1,decalreflet);
    reflet = [decalageReflet,RefRaccourci]*attenuation + Ref;
    FFTref = fftshift(fft(Ref));
    FFTcopie = fftshift(fft(reflet));
    vec = fftshift(ifft(conj(FFTref).*FFTcopie));
end

% function [reflet] = bulbizarre(Ref,decalreflet,attenuation)
%     L = length(Ref);
%     removeEnd = L - decalreflet+1;
%     RefRaccourci = Ref;
%     RefRaccourci(removeEnd:end) = [];
%     decalageReflet = zeros(1,decalreflet);
%     reflet = [decalageReflet,RefRaccourci]*attenuation;
% 
% end


function [vec] = fourier(r)
    L = length(r);
    R = fftshift(fft(r));
    vec = abs(R/L);
end

 

function [vec] = fourier_inverse(R)
    %L = length(R);
    r = fftshift(ifft(R));
    vec = abs(r);
end
