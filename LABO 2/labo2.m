clear all; close all; clc;

% Pour enlever le message d'erreur.
MSGID = 'signal:hilbert:Ignore';
warning('off', MSGID);

load("data_labo_reflexion.mat");

% Donne la plage du signal qui va être analysée
global precision sizeComparaison;
precision = 110000;
sizeComparaison = 1500;

global ref REF;
ref = Reference;
REF = fftshift(fft(ref));
[~,ind_max_corr] = max(Corr7);
corr = abs(hilbert(Corr7(ind_max_corr-sizeComparaison/2:ind_max_corr+sizeComparaison/2-1)));
[delay, attenuation,compare_fictif,compare_corr] = findParameters(corr);
delay = 87;
attenuation = 0.71;
%plotGraphs(compare_corr, compare_fictif);
%vec = signalSynth(delay, attenuation);



function [delay,attenuation,compare_fictif_end,compare_corr_end] = findParameters(corr)
    %global sizeComparaison;

    delays = 1:100;
    attenuations = 0.3:0.1:0.8;
    D = repmat(delays', length(attenuations),1);
    A = repelem(attenuations',length(delays));

    errs = zeros(1,length(delays)*length(attenuations));
    fictif = signalSynth(D, A);
    
    [max_fict,~] = max(fictif,[],2);
    [max_corr,~] = max(corr);

    rapport = max_corr/max_fict;
    compare_fictif = fictif.*rapport';
    
    
    % La méthode de comparaison de l'efficacité va se baser sur la MSE.
    for i = 1:length(errs)
        errs(i) = immse(compare_fictif(i,:), corr(1,:));
    end

    %mesh(attenuations, delays, errs);
    [~, index_min_errs] = min(errs(:));
    %[index_del, index_att] = ind2sub(size(errs), index_min_errs);

    delay = D(index_min_errs);
    attenuation = A(index_min_errs);
    compare_fictif_end = compare_fictif(index_min_errs,:);
    compare_corr_end = corr(1,:);
end


% Cree le signal synthetique avec :
% D : le delay 
% A : l'attenuation
function [vec] = signalSynth(D, A)
    global ref REF sizeComparaison;
    indices = length(D);
    
    a = zeros(indices, length(ref));
    for i = 1:length(D)
        a(i,D(i):end) = ref(1,1:end-D(i)+1)*A(i);
    end
    
    %{
    reflet = repmat(ref, indices,1) + a;
    REFLET = fftshift(fft(reflet));
    vec = (ifft(ifftshift(conj(REF).*REFLET)));
    [~, ind_max] = max(vec,[],2);
    vec = vec(:,ind_max-sizeComparaison/2:ind_max+sizeComparaison/2-1);
    %}
    reflet = repmat(ref, indices,1) + a;
    
    e = fft(reflet);
    f = fft(ref);
    vec = ifftshift(ifft(e.*conj(f)));
    
    
    REFLET = fft(reflet,[],1);
    
    figure()
    plot(abs(vec(20,:)));

    
    %figure()
    %plot(abs(REF))
    
    vec = (ifft(ifftshift(conj(REF).*REFLET)));
    [~, ind_max] = max(vec,[],2);
    vec = vec(:,ind_max-sizeComparaison/2:ind_max+sizeComparaison/2-1);
    
    %figure()
    %plot(vec2(20,:))
    %hold on;
    %plot(abs(hilbert(vec(20,:))),'--');
    
    vec = abs(hilbert(vec));
end


function [] = plotGraphs(compare_corr, compare_fictif)

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


