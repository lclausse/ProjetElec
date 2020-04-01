clear;
clc;

load("data_labo_reflexion.mat");

subplot(2, 2, 1);
plot(abs(Corr7));
decalagedebase = 102250;
decalsupp = 115;
attenuation = 0.6666;
test1 = miaouss(decalagedebase,decalsupp,attenuation,Reference);
subplot(2,2,2);
plot(abs(test1));
test2 = tortank(decalagedebase+1100,attenuation,Reference);
subplot(2,2,3);
plot(abs(test2));
test3 = tortank(decalagedebase+1200,attenuation,Reference);
subplot(2,2,4);
plot(abs(test3));






function [jeanphilippe] = tortank(sampledecal,attenuation,Ref)
    decalage = zeros(1,sampledecal);
    refcopie = [decalage,Ref]*attenuation;
    nouvRef = [Ref,decalage];
    FFTref = fourier(nouvRef);
    FFTcopie = fourier(refcopie);
    jeanphilippe = fourier_inverse(conj(FFTref).*FFTcopie);
end


function [vec] = miaouss(decalbase,decalreflet,attenuation,Ref)
    decalage = zeros(1,decalbase);
    decalageReflet = zeros(1,decalreflet);
    nouvRef = [decalage,Ref,decalageReflet];
    refCopie = [decalage,decalageReflet,Ref]*attenuation;
    FFTref = fourier(nouvRef);
    FFTcopie = fourier(refCopie);
    vec = fourier_inverse(conj(FFTref).*FFTcopie);


end


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
