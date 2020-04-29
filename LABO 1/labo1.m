clear all;
close all;
clc;
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
TDOARef = Time.PcAtRoutineRun(1:end);    %Timestamps of the PC clock when it enters its localisation routine
c = Time.PcAtTsFrame; c = flip(c); %Timestamps of the PC clock when it receive a measurement from the total station
d = Time.TsAtTsFrame; d = flip(d); %Timestamps of the TotalStation clock when it takes the measurement
PcAtRoutineRunVect = seconds(TDOARef-TDOARef(end));
PcAtTsFrameVect = seconds(c-TDOARef(end)); 
TsAtTsFrameVect = seconds(d-TDOARef(end));
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
global xRef;
xRef = xCalTag(1,1);
global yRef;
yRef = xCalTag(2,1);


%scatter3(xTotalStationSync(1,:), xTotalStationSync(2,:), xTotalStationSync(3,:), 'filled');
%plot(RawSignalRx1(1,:)-RawSignalRx1(2,:))

% Tableau de 20x81920 -> 20 points avec 81920 samples pour chaque et
% une fréquence d'échantillonage de FsReference

i = 1;
% On enlève la moyenne des signaux (partie DC)
r1 = RawSignalRx1(i,:) - mean(RawSignalRx1(i,:));
r2 = RawSignalRx2(i,:) - mean(RawSignalRx2(i,:));
r3 = RawSignalRx3(i,:) - mean(RawSignalRx3(i,:));
r4 = RawSignalRx4(i,:) - mean(RawSignalRx4(i,:));

% Limites où on coupe les signaux pour isoler partie de la balise et ref
limiteTemp = 31000;
limiteLowRef = 36000;
limiteHighRef = 80000;

% On coupe les signaux
r1Balise = r1(1:limiteTemp);
r1Ref = r1(limiteLowRef:limiteHighRef);
r2Balise = r2(1:limiteTemp);
r2Ref = r2(limiteLowRef:limiteHighRef);
r3Balise = r3(1:limiteTemp);
r3Ref = r3(limiteLowRef:limiteHighRef);
r4Balise = r4(1:limiteTemp);
r4Ref = r4(limiteLowRef:limiteHighRef);

global Fs Fsample tau;
% Fréquence d'échantillonnage
Fsample = FsReference;
% Fréquence du signal
Fs = FsRawSignal;

% Les TDOA qu'on devrait trouver. Calcul géométrique.
TDOARefExact = trueTDOARef();


%{
[fRf,R1Balise] = toRf(r1Balise,Fs);
[fRf,R2Balise] = toRf(r2Balise,Fs);
[fRf,R3Balise] = toRf(r3Balise,Fs);
[fRf,R4Balise] = toRf(r4Balise,Fs);

[fRef,R1Ref] = toRf(r1Ref,Fs);
[fRef,R2Ref] = toRf(r2Ref,Fs);
[fRef,R3Ref] = toRf(r3Ref,Fs);
[fRef,R4Ref] = toRf(r4Ref,Fs);

R1 = abs(fftshift(fft(r1))/length(r1));
R2 = abs(fftshift(fft(r2))/length(r2));
R3 = abs(fftshift(fft(r3))/length(r3));
R4 = abs(fftshift(fft(r4))/length(r4));
%}
%{
figure()
subplot(2,1,1)
plot(r1)
subplot(2,1,2)
plot(r2)


figure();
subplot(3,2,1);
plot(r1(1:limiteTemp));
subplot(3,2,2);
plot(r2(1:limiteTemp));

subplot(3,2,3);
plot(fRf,R1Balise);
subplot(3,2,4);
plot(fRf,R2Balise);

subplot(3,2,5);
%}
%{
figure();
subplot(2,1,1);
corr = fourier_inverse(conj(R1) .* R2);
plot(corr);
%}

%xEst = zeros(1,length(TDOA));
%yEst = zeros(1,length(TDOA));
%resnorm = zeros(1,length(TDOA));




% Delais des balises
TDOABal = [findDelay(r1Balise,r2Balise),...
           findDelay(r1Balise,r3Balise),...
           findDelay(r1Balise,r4Balise),...
           findDelay(r3Balise,r2Balise),...
           findDelay(r4Balise,r2Balise),...
           findDelay(r3Balise,r4Balise)];
        
% Delais de référence
TDOARef = [findDelay(r1Ref,r2Ref),...
           findDelay(r1Ref,r3Ref),...
           findDelay(r1Ref,r4Ref),...
           findDelay(r3Ref,r2Ref),...
           findDelay(r4Ref,r2Ref),...
           findDelay(r3Ref,r4Ref)];
        
% Correction des erreurs statiques
errCorrection = TDOARef - TDOARefExact;

% !!!!!!!!!!!!!!!! IL Y A DU CACA DE SIGNES ICI !!!!!!!!!!!!!!!!

% Delais corrigés des balises
tau = TDOABal - errCorrection;

x0 = [5,5];
[temp,resnorm] = lsqnonlin(@func,x0);
xEst = temp(1);
yEst = temp(2);


figure()
scatter(1:length(tau),tau,'filled');
hold on;
scatter(1:length(TDOARef),TDOARef);
hold on;
scatter(1:length(TDOARefExact),TDOARefExact);
grid on;
legend('Corrected TDOA','Computed TDOA','Geometric TDOA');


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
xlabel('position [m]');
ylabel('position [m]');
grid on;
%}



% Fonction qui donne le TDOA entre les deux signaux.
% Calcule la corrélation puis le TDOA.
function timeDelay = findDelay(r1,r2)
    global Fsample;
    
    % On décale les signaux pour centrer la corrélation
    if (length(r1) < length(r2))
        c = [ zeros(1,length(r2)-1) r1 zeros(1,length(r2)-length(r1)) ];
        d = [ r2 zeros(1,length(r2)-1) ];
    else
        c = [ zeros(1,length(r1)-1) r1 ];
        d = [ r2 zeros(1,length(r1)-length(r2)+length(r1)-1) ];
    end
    
    e = fft(c);
    f = fft(d);
    g = e.*conj(f);
    h = ifft(g);
    
    % Axe x de la corrélation centrée
    samplesScale = -length(h)/2:1:length(h)/2 - 1;
    % Maximum de cette corrélation -> index
    [~,maxIndex] = max(abs(h));
    % On prend cet index sur l'axe x trouvé avant
    sampleDelay = samplesScale(maxIndex);
    % On calcule un temps avec la fréquence
    timeDelay = sampleDelay/Fsample;
    
    % Plot pour le rapport
    %{
    figure()
    subplot(4,1,1);
    plot(r1);
    title('Rx1');
    
    subplot(4,1,2);
    plot(r2);
    title('Rx2');
    
    subplot(4,1,3);
    plot(c);
    title('Rx1 décalé');
    
    subplot(4,1,4);
    plot(d);
    title('Rx2 décalé');
    xlabel('samples');
    
    figure()
    plot(samplesScale,h)
    title('Correlation');
    xlabel('samples');
    %}
    
end


% TDOA qu'on devrait avoir pour la balise de référence
% On part des positions. Constant pour tous les points.
function TDOARef = trueTDOARef()
    global lightSpeed x1 x2 y1 y2 xRef yRef;
    tau = zeros(1,6);
    tau(1) = (sqrt((xRef-x2(1))^2+(yRef-y2(1))^2)-sqrt((xRef-x1(1))^2+(yRef-y1(1))^2)) / lightSpeed;
    tau(2) = (sqrt((xRef-x2(2))^2+(yRef-y2(2))^2)-sqrt((xRef-x1(2))^2+(yRef-y1(2))^2)) / lightSpeed;
    tau(3) = (sqrt((xRef-x2(3))^2+(yRef-y2(3))^2)-sqrt((xRef-x1(3))^2+(yRef-y1(3))^2)) / lightSpeed;
    tau(4) = (sqrt((xRef-x2(4))^2+(yRef-y2(4))^2)-sqrt((xRef-x1(4))^2+(yRef-y1(4))^2)) / lightSpeed;
    tau(5) = (sqrt((xRef-x2(5))^2+(yRef-y2(5))^2)-sqrt((xRef-x1(5))^2+(yRef-y1(5))^2)) / lightSpeed;
    tau(6) = (sqrt((xRef-x2(6))^2+(yRef-y2(6))^2)-sqrt((xRef-x1(6))^2+(yRef-y1(6))^2)) / lightSpeed;
    TDOARef = tau;
end


% Fonction non linéaire à résoudre pour trouver les positions
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






% ------------------------------------------------------------------------
% --------------------------- bullshit -----------------------------------
% ------------------------------------------------------------------------

 
% Fonction du net
%{
function t = findCorrAndDelay(r1,r2)
    global Fs;
    [acor,lag] = xcorr(r2,r1);
    [~,I] = max(abs(acor));
    lagDiff = lag(I);
    t = lagDiff/Fs;
    %a3 = gca;
    %a3.XTick = sort([-3000:1000:3000 lagDiff])
end
%}  

%{
function [vec] = fourier(r)
    L = length(r);
    R = fftshift(fft(r));
    vec = abs(R/L);
end

function [vec] = fourier_inverse(R)
    r = fftshift(ifft(R));
    vec = real(r);
end
%}

%{
% Passe le signal vers de la radio fréquence.
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
%}
