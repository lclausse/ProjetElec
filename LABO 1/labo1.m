clear
clc
%importfile('Data_Measured.mat')
%importfile('Data_Synthetic.mat')
load('Data_Lab1_2.mat')


global lightSpeed;
lightSpeed = 299792458;

% ------------------------------------------------------------------------
% ---------------------------- CODE DES PROFS ----------------------------
% ------------------------------------------------------------------------
%The piece of code below Interpolate the Total Station measurements at the
%time the UWB signals are acquired.
%The input variable (Time and xTotalStation) are supposed to be loaded from
%the files saved during the lab
%Due to the way the variable are updated on the main computer during the
%lab, some variable have to be flipped.
%We consider that the time of the first localisation run on the PC (==a(end)) is t=0.
a = Time.PcAtRoutineRun(1:end);    %Timestamps of the PC clock when it enters its localisation routine
c = Time.PcAtTsFrame; c = flip(c); %Timestamps of the PC clock when it receive a measurement from the total station
d = Time.TsAtTsFrame; d = flip(d); %Timestamps of the TotalStation clock when it takes the measurement
PcAtRoutineRunVect = seconds(a-a(end));
PcAtTsFrameVect = seconds(c-a(end)); 
TsAtTsFrameVect = seconds(d-a(end));
ClockOffset = seconds(d(end)-c(end)); %Offset between the PC clock and the Total Station clock
TsAtTsFrameVect = TsAtTsFrameVect - ClockOffset; 
xTsFlip = flip(xTotalStation,2);
xTotalStationSync = spline(TsAtTsFrameVect,xTsFlip,PcAtRoutineRunVect+0.2035); %The constant delay 0.2035s is based on error minimisation on a large data set.
% ------------------------------------------------------------------------
% --------------------------------- FIN ----------------------------------
% ------------------------------------------------------------------------

global x1;
x1 = [xReceivers(1,1), xReceivers(1,1), xReceivers(1,1), xReceivers(1,2), xReceivers(1,2),xReceivers(1,3)];
global y1;
y1 = [xReceivers(2,1), xReceivers(2,1), xReceivers(2,1), xReceivers(2,2), xReceivers(2,2),xReceivers(2,3)];
global x2;
x2 = [xReceivers(1,2), xReceivers(1,3), xReceivers(1,4), xReceivers(1,3), xReceivers(1,4),xReceivers(1,4)];
global y2;
y2 = [xReceivers(2,2), xReceivers(2,3), xReceivers(2,4), xReceivers(2,3), xReceivers(2,4),xReceivers(2,4)];

%scatter3(xTotalStationSync(1,:), xTotalStationSync(2,:), xTotalStationSync(3,:), 'filled');
%plot(RawSignalRx1(1,:)-RawSignalRx1(2,:))

% Tableau de 20x81920 -> 20 points avec 81920 samples pour chaque et
% une fréquence d'échantillonage de FsReference

r1 = RawSignalRx1(1,:) - mean(RawSignalRx1(1,:));
r2 = RawSignalRx2(1,:) - mean(RawSignalRx2(1,:));
r3 = RawSignalRx3(1,:) - mean(RawSignalRx3(1,:));
r4 = RawSignalRx4(1,:) - mean(RawSignalRx4(1,:));

% fréquence d'échantillonage :
Fsample = FsReference;
Fs = FsRawSignal;
%L = length(r1);
%f = (-L/2:L/2-1)*(Fsample/L);

limiteTemp = 30000;

[f,R1] = toRf(r1(1:limiteTemp),Fs);
[f,R2] = toRf(r2(1:limiteTemp),Fs);
[f,R3] = toRf(r3(1:limiteTemp),Fs);
[f,R4] = toRf(r4(1:limiteTemp),Fs);

figure();
subplot(3,2,1);
plot(r1(1:limiteTemp));
subplot(3,2,2);
plot(r2(1:limiteTemp));

subplot(3,2,3);
plot(f,R1);
subplot(3,2,4);
plot(f,R2);

subplot(3,2,5);
corr = fourier_inverse(conj(R1) .* R2);
plot(corr);


%{
R1 = fourier(r1);
R2 = fourier(r2);
R3 = fourier(r3);
R4 = fourier(r4);

figure();
subplot(3,2,1);
plot(r1)
title('Signal r_1')
xlabel('samples')
ylabel('r_1')

subplot(3,2,2);
plot(r2)
title('Signal r_2')
xlabel('samples')
ylabel('r_2')

subplot(3,2,3);
plot(f, R1)
title('Spectre de r_1')
xlabel('f [Hz]')
ylabel('|R1|')

subplot(3,2,4);
plot(f, R2)
title('Spectre de r_2')
xlabel('f [Hz]')
ylabel('|R2|')

subplot(3,2,5);
CORR = conj(R1) .* R2;
corr = fourier_inverse(CORR);
plot(corr)
title('Corrélation r_1 r_2')

[val,ind] = max(corr)
ind / FsReference
%}

% ---- PLOT ----
%{
figure();
scatter(xReceivers(1,:),xReceivers(2,:),'filled') % Recepteurs
str = [" R1"," R2"," R3"," R4"];
text(xReceivers(1,:),xReceivers(2,:), str)      
hold on;
scatter(xCalTag(1,1),xCalTag(2,1), 'filled'); % Balise de calibrage
text(xCalTag(1,1),xCalTag(2,1), ' Ref')
hold on;
%scatter(xEst,yEst,'*')                          % Estimation de la position
hold on;
%text(-8, 12, "MSE = " + err)
hold on;
scatter(xTotalStationSync(1,:),xTotalStationSync(2,:),'k','.') % Position au laser
legend('Récepteurs','Référence','True pos','Location','north')
%}


function F = func(x)
    global x1 x2 y1 y2 lightSpeed tau;
    f1 = sqrt((x(1)-x2(1))^2+(x(2)-y2(1))^2)-sqrt((x(1)-x1(1))^2+(x(2)-y1(1))^2)-lightSpeed*tau(1);
    f2 = sqrt((x(1)-x2(2))^2+(x(2)-y2(2))^2)-sqrt((x(1)-x1(2))^2+(x(2)-y1(2))^2)-lightSpeed*tau(2);
    f3 = sqrt((x(1)-x2(3))^2+(x(2)-y2(3))^2)-sqrt((x(1)-x1(3))^2+(x(2)-y1(3))^2)-lightSpeed*tau(3);
    f4 = sqrt((x(1)-x2(4))^2+(x(2)-y2(4))^2)-sqrt((x(1)-x1(4))^2+(x(2)-y1(4))^2)-lightSpeed*tau(4);
    f5 = sqrt((x(1)-x2(5))^2+(x(2)-y2(5))^2)-sqrt((x(1)-x1(5))^2+(x(2)-y1(5))^2)-lightSpeed*tau(5);
    f6 = sqrt((x(1)-x2(6))^2+(x(2)-y2(6))^2)-sqrt((x(1)-x1(6))^2+(x(2)-y1(6))^2)-lightSpeed*tau(6);
    F = [f1, f2, f3, f4, f5, f6];
end

function [vec] = fourier(r)
    L = length(r);
    R = fftshift(fft(r));
    vec = abs(R/L);
end

function [vec] = fourier_inverse(R)
    r = fftshift(ifft(R));
    vec = real(r);
end

% r  : signal temporel à mettre en radio-fréquences
% Fs : Fréquence du signal
% Fsample : Fréquence d'échantillonage
function [f,R] = toRf(r, Fs)
    numZeros = 3;
    L_per = numZeros * length(r) - (numZeros-1);
    r_per = zeros(1,L_per);
    r_per(1:numZeros:end) = r;
    R_per = abs(fftshift(fft(r_per))/L_per);
    f_per = numZeros*(-L_per/2:L_per/2-1)*(Fs/L_per);
    
    [val,indexFs] = min(abs(f_per-Fs));
    R_per(L_per-indexFs:indexFs) = 0;
    f = f_per;
    R = R_per;
end
