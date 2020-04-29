clear
clc

load('Data_pour_notre_groupe.mat')

Fs = FsRawSignal;

r1 = RawSignalRx1 - mean(RawSignalRx1); % Enlever moyenne DC
r2 = RawSignalRx2 - mean(RawSignalRx2);
L = length(r1);
if L ~= length(r2)
    disp("Error : the two signals have different size.");
end


Tsample = 1 / FsReference;
time = Tsample * (0:L-1);

subplot(4,2,1);
plot(time, r1);
title('Signal temporel Rx1')
xlabel('t [s]')
ylabel('r_1(t)')

subplot(4,2,2);
plot(time, r2);
title('Signal temporel Rx2')
xlabel('t [s]')
ylabel('r_2(t)')


R1 = fourier(r1);
R2 = fourier(r2);

f = (-L/2:L/2-1)*(Fs/L); % zero-centered frequency range

subplot(4,2,3);
plot(f,R1) 
title('Spectre de r_1')
xlabel('f [Hz]')
ylabel('|R_1(f)|')

subplot(4,2,4);
plot(f,R2) 
title('Spectre de r_2')
xlabel('f [Hz]')
ylabel('|R_2(f)|')


L_per = 3 * length(r1) - 2;% * (length(r1) - 1);

f_per = 3*(-L_per/2:L_per/2-1)*(Fs/L_per);
% On ajoute des 0 entre les samples pour périodiser le spectre
r1_per = zeros(1,L_per);
r1_per(1:3:end) = r1;
R1_per = fourier(r1_per);

r2_per = zeros(1,L_per);
r2_per(1:3:end) = r2;
R2_per = fourier(r2_per);

subplot(4,2,5);
plot(f_per,R1_per)
title('Spectre périodisé de r_1')
xlabel('f [Hz]')
ylabel('|R1_{per}(f)|')

subplot(4,2,6);
plot(f_per,R2_per) 
title('Spectre périodisé de r_2')
xlabel('f [Hz]')
ylabel('|R2_{per}(f)|')

% On va mtn filtrer ce qui est en dessous de Fs = 3.2 GHz
posFs = cast((2/3) * L_per, 'int32');
f_per(posFs)

R1_per_filtre = R1_per;
R1_per_filtre(1:posFs) = 0;

R2_per_filtre = R2_per;
R2_per_filtre(1:posFs) = 0;

subplot(4,2,7);
plot(f_per,R1_per_filtre)
title('Spectre périodisé filtré de r_1')
xlabel('f [Hz]')
ylabel('|R1_{filt}(f)|')

subplot(4,2,8);
plot(f_per,R2_per_filtre) 
title('Spectre périodisé filtré de r_2')
xlabel('f [Hz]')
ylabel('|R2_{filt}(f)|')


% Corrélation entre signal r1 et r2 : 
sample_corr = (-L_per/2:L_per/2-1); %(Tsample/3) * ;

% En mode RF -> Remis bonne fréquence
CORR_RF = conj(R1_per_filtre) .* R2_per_filtre;
corr_RF = fourier_inverse(CORR_RF);

% En mode bande de base
CORR_BB = conj(R1) .* R2;
corr_BB = fourier_inverse(CORR_BB);

figure();

subplot(2,1,1);
plot(CORR_RF);

subplot(2,1,2);
plot(sample_corr, corr_RF);


function [vec] = fourier(r)
    L = length(r);
    R = fftshift(fft(r));
    vec = abs(R/L);
end

function [vec] = fourier_inverse(R)
    %L = length(R);
    r = fftshift(ifft(R));
    vec = real(r);
end





















