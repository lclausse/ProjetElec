clear all; close all; clc;
% Pour enlever le message d'erreur.
MSGID = 'signal:hilbert:Ignore';
warning('off', MSGID);

load("data_labo_reflexion.mat");

global sizeComparaison ref REF;
sizeComparaison = 500;
ref = Reference;
REF = fftshift(fft(ref));

corr = abs(hilbert(Corr7));
[max_corr, i_max_corr] = max(corr);
corr_compare = corr(i_max_corr-sizeComparaison/2:i_max_corr+sizeComparaison/2-1);

[delay, attenuation, fictif] = optimisation(corr_compare);

% --- PLOT RAPPORT ---
t = delay/(38.4 * 10^9);
d = t * 3 * 10^8;
figure()
plot(fictif);
hold on;
plot(corr_compare);
grid on;
legend('Signal fictif','Signal corrélé');
title('Comparaison signal fictif et corrélé');
% --- END PLOT RAPPORT ---

function [delay, attenuation, fictif] = optimisation(corr)
    global sizeComparaison;
    delays = 50:100;
    attenuations = 0.1:0.01:0.9;
    D = repmat(delays', length(attenuations),1);
    A = repelem(attenuations',length(delays));
    
    errs = zeros(1,length(delays)*length(attenuations));
    max_corr = max(corr);
    
    for i = 1:length(errs)
        fict = signalSynth(D(i), A(i));
        [max_fict, i_max_fict] = max(fict);  
        rapport = max_fict/max_corr;
        compare_fictif = fict(i_max_fict-sizeComparaison/2:i_max_fict+sizeComparaison/2-1)/rapport;
        errs(i) = immse(compare_fictif, corr);
    end
    
    [~, index_min_errs] = min(errs);
    delay       = D(index_min_errs);
    attenuation = A(index_min_errs);
    
    fictif = signalSynth(delay, attenuation);
    [max_fict, i_max_fict] = max(fictif);
    rapport = max_fict/max_corr;
    fictif = fictif(i_max_fict-sizeComparaison/2:i_max_fict+sizeComparaison/2-1)/rapport;
end

function [vec] = signalSynth(D,A)
    global ref REF;
    a = zeros(1, length(ref));
    a(1,D:end) = ref(1,1:end-D+1)*A;
    reflet = ref + a;
    REFLET = fftshift(fft(reflet));
    vec = abs(hilbert(ifftshift(ifft(conj(REF).*REFLET))));
end


