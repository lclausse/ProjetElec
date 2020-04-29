clear;
clc;

% Pour enlever le message d'erreur.
MSGID = 'signal:hilbert:Ignore';
warning('off', MSGID);

load("data_labo_reflexion.mat");

% Donne la plage du signal qui va être analysée
global precision;
precision = 10000;
global sizeComparaison;
sizeComparaison = 1500;


corr = abs(hilbert(Corr7(1:precision)));
ref = Reference(1:precision);

REF = fftshift(fft(ref));
[delay, attenuation] = findParameters(corr, ref, REF);
%delay = 87;
%attenuation = 0.71;
%plotGraphsRapport(corr, Reference, delay, attenuation);


function [delay, attenuation] = findParameters(corr, ref, REF)
    global sizeComparaison;
    global precision;

    delays = 70:100;
    attenuations = 0.6:0.01:0.8;
    D = repmat(delays', length(attenuations), 1);
    A = repelem(attenuations', length(delays));
    
    errs = zeros(length(delays)*length(attenuations));
    fictif = signalSynth(D, A, ref, REF);
    

    % Il faut alligner les deux vecteurs sur la pique.
    [max_fictif, indice_max_fictif] = max(fictif,[],2);
    [max_corr, indice_max_corr] = max(corr);

    rapport = max_corr/max_fictif;
    
    compare_fictif = fictif(:,indice_max_fictif-sizeComparaison/2:indice_max_fictif+sizeComparaison/2-1);%.*rapport';
    compare_corr = corr(indice_max_corr-sizeComparaison/2:indice_max_corr+sizeComparaison/2-1);
    
    indice_max_fictif(200)
    subplot(2,1,1);
    plot(compare_fictif(200,:))
    subplot(2,1,2);
    plot(compare_corr)
    

    % La méthode de comparaison de l'efficacité va se baser sur la MSE.
    % Modèle de régression linéaire à deux variables (delay, attenuation).

    errs(index_del,index_att) = immse(compare_fictif, compare_corr);


    %mesh(attenuations, delays, errs);
    [min_errs, index_min_errs] = min(errs(:));
    [index_del, index_att] = ind2sub(size(errs), index_min_errs);

    delay = delays(index_del);
    attenuation = attenuations(index_att);

end



function [vec] = signalSynth(D, A, ref, REF)
    global precision;
    indices = length(D);
    
    a = zeros(indices, precision);
    for i = 1:length(D) % dommage que pas vectorisé, mais matrice pas rectangulaire...
        a(i,D(i):end) = ref(1,1:end-D(i)+1)*A(i);
    end
        
    reflet = repmat(ref, indices,1) + a;    
    REFLET = fftshift(fft(reflet));
    vec = abs(hilbert(fftshift(ifft(conj(REF).*REFLET))));
end








function [] = plotGraphs(correlation, ref, delay, attenuation)
    corr = abs(hilbert(correlation));
    fictif = signalSynth(delay, attenuation, ref);
    global sizeComparaison;

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


function [] = plotGraphsRapport(correlation, ref, delay, attenuation)
    corr = abs(hilbert(correlation));
    fictif = signalSynth(delay, attenuation, ref);
    global sizeComparaison;

    [max_fictif, indice_max_fictif] = max(fictif);
    [max_corr, indice_max_corr] = max(corr);

    rapport = max_fictif/max_corr;

    figure();

    subplot(3,2,1);
    plot(fictif);
    title('Signal fictif')

    subplot(3,2,2);
    plot(corr);
    title('Signal corrélé')

    compare_fictif = fictif(indice_max_fictif-sizeComparaison/2:indice_max_fictif+sizeComparaison/2-1);
    compare_corr = corr(indice_max_corr-sizeComparaison/2:indice_max_corr+sizeComparaison/2-1)*rapport; % pour le moment je les match comme ça

    subplot(3,2,3);
    plot(compare_fictif);
    title('Signal fictif centré')

    subplot(3,2,4);
    plot(compare_corr);
    title('Signal corrélé centré')

end


