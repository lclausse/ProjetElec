clear;
clc;

load("data_labo_reflexion.mat");

[delay, attenuation] = optimisation(Corr7, Reference)
plotGraphs(Corr7, Reference, delay, attenuation);


function [delay, attenuation] = optimisation(correlation, ref)
    delays = 70:100;
    attenuations = 0.6:0.001:0.8;
    errs = zeros(length(delays),length(attenuations));
    
    corr = abs(hilbert(correlation));
    sizeComparaison = 500;
        
    %{ 
    % A vectoriser!!
    for index_del = 1:length(delays)
        for index_att = 1:length(attenuations)
            
            fictif = racaillou(delays(index_del), attenuations(index_att), ref);
            % Il faut alligner les deux vecteurs sur la pique.
            [max_fictif, indice_max_fictif] = max(fictif);
            [max_corr, indice_max_corr] = max(corr);
    
            rapport = max_fictif/max_corr;
    
            compare_fictif = fictif(indice_max_fictif-sizeComparaison/2:indice_max_fictif+sizeComparaison/2-1);
            compare_corr = corr(indice_max_corr-sizeComparaison/2:indice_max_corr+sizeComparaison/2-1)*rapport; % pour le moment je les match comme ça
    
            % La méthode de comparaison de l'efficacité va se baser sur la MSE.
            % Modèle de régression linéaire à deux variables (delay, attenuation). 
    
            errs(index_del,index_att) = immse(compare_fictif, compare_corr);
        	
        end
    end
    %}
    
    for index_del = 1:length(delays)
        for index_att = 1:length(attenuations)
            
            fictif = racaillou(delays(index_del), attenuations(index_att), ref);
            % Il faut alligner les deux vecteurs sur la pique.
            [max_fictif, indice_max_fictif] = max(fictif);
            [max_corr, indice_max_corr] = max(corr);
    
            rapport = max_fictif/max_corr;
    
            compare_fictif = fictif(indice_max_fictif-sizeComparaison/2:indice_max_fictif+sizeComparaison/2-1);
            compare_corr = corr(indice_max_corr-sizeComparaison/2:indice_max_corr+sizeComparaison/2-1)*rapport; % pour le moment je les match comme ça
    
            % La méthode de comparaison de l'efficacité va se baser sur la MSE.
            % Modèle de régression linéaire à deux variables (delay, attenuation). 
    
            errs(index_del,index_att) = immse(compare_fictif, compare_corr);
        	
        end
    end
    
    
    %mesh(attenuations, delays, errs);
    [min_errs, index_min_errs] = min(errs(:));
    [index_del, index_att] = ind2sub(size(errs), index_min_errs);
    
    delay = delays(index_del);
    attenuation = attenuations(index_att);

end


function [] = plotGraphs(correlation, ref, delay, attenuation)
    corr = abs(hilbert(correlation));
    fictif = racaillou(delay, attenuation, ref);
    sizeComparaison = 500;
    
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




function [vec] = racaillou(decalreflet,attenuation,Ref)
    L = length(Ref);
    removeEnd = L - decalreflet+1;
    RefRaccourci = Ref;
    RefRaccourci(removeEnd:end) = [];
    decalageReflet = zeros(1,decalreflet);
    reflet = [decalageReflet,RefRaccourci]*attenuation + Ref;
    FFTref = fftshift(fft(Ref));
    FFTcopie = fftshift(fft(reflet));
    vec = abs(hilbert(fftshift(ifft(conj(FFTref).*FFTcopie))));
end
