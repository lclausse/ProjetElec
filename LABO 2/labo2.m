clear all; close all; clc;

% Pour enlever le message d'erreur.
MSGID = 'signal:hilbert:Ignore';
warning('off', MSGID);

load("data_labo_reflexion.mat");

global sizeComparaison;
sizeComparaison = 1500;

global ref REF;
ref = Reference;
REF = fftshift(fft(ref));


%subplot(2, 2, 1);
%plot(abs(hilbert(Corr7)));
%decalsupp = 85;
%attenuation = 0.66;


[delay, attenuation] = optimisation(Corr7);
%plotGraphs(Corr7, Reference, delay, attenuation);

%{
test1 = signalSynth(113,attenuation,Reference);
subplot(2,2,2);
plot((test1));
test2 = signalSynth(115,attenuation,Reference);
subplot(2,2,3);
plot((test2));
test3 = signalSynth(118,attenuation,Reference); %118 ca donne bien
subplot(2,2,4);
plot((test3));
%}

function [delay, attenuation] = optimisation(correlation)
    delays = 50:100;
    attenuations = 0.1:0.1:0.9;
    D = repmat(delays', length(attenuations),1);
    A = repelem(attenuations',length(delays));
    
    errs = zeros(1,length(delays)*length(attenuations));
    
    corr = abs(hilbert(correlation));
    sizeComparaison = 1500;
        
    for i = 1:length(errs)
        
        fict = signalSynth(D(i), A(i));
        % Il faut alligner les deux vecteurs sur la pique.
        [max_fict, i_max_fict] = max(fict);
        [max_corr, i_max_corr] = max(corr);

        rapport = max_fict/max_corr;

        compare_fictif = fict(i_max_fict-sizeComparaison/2:i_max_fict+sizeComparaison/2-1);
        compare_corr = corr(i_max_corr-sizeComparaison/2:i_max_corr+sizeComparaison/2-1)*rapport;

        errs(i) = immse(compare_fictif, compare_corr);
    end
    
    
    %mesh(attenuations, delays, errs);
    [~, index_min_errs] = min(errs);
    delay = D(index_min_errs);
    attenuation = A(index_min_errs);
    
    fictifEnd = signalSynth(delay, attenuation);
    [max_fict, i_max_fict] = max(fictifEnd);
    [max_corr, i_max_corr] = max(corr);
    rapport = max_fict/max_corr;
    figure();
    plot(fictifEnd(i_max_fict-sizeComparaison/2:i_max_fict+sizeComparaison/2-1));
    hold on;
    plot(corr(i_max_corr-sizeComparaison/2:i_max_corr+sizeComparaison/2-1)*rapport);

end


function [] = plotGraphs(correlation, ref, delay, attenuation)
    corr = abs(hilbert(correlation));
    fictif = signalSynth(delay, attenuation);
    sizeComparaison = 1000;
    
    [max_fictif, indice_max_fictif] = max(fictif);
    [max_corr, indice_max_corr] = max(corr);
    
    rapport = max_fictif/max_corr;
    
    compare_fictif = fictif(indice_max_fictif-sizeComparaison/2:indice_max_fictif+sizeComparaison/2-1);
    compare_corr = corr(indice_max_corr-sizeComparaison/2:indice_max_corr+sizeComparaison/2-1)*rapport; % pour le moment je les match comme ça
    
    figure();
    subplot(3,1,1);
    plot(compare_corr);
    
    subplot(3,1,2);
    plot(compare_fictif);
    
    subplot(3,1,3);
    plot(compare_corr - compare_fictif);
end


function [vec] = signalSynth(D,A)
    global ref REF;
    a = zeros(1, length(ref));
    a(1,D:end) = ref(1,1:end-D+1)*A;
    reflet = ref + a;
    REFLET = fftshift(fft(reflet));
    vec = abs(hilbert(ifftshift(ifft(conj(REF).*REFLET))));
end


