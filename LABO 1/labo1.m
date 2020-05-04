clear all;
close all;
clc;
%importfile('Data_Measured.mat')
%importfile('Data_Synthetic.mat')
load('Data_Lab1_3.mat')

% Pour enlever le message d'erreur.
MSGID = 'MATLAB:declareGlobalBeforeUse';
warning('off', MSGID);
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

global Fs Fsample FsRawSignal tau;
% Fréquence d'échantillonnage
Fsample = FsReference;
% Fréquence du signal
Fs = 3*FsRawSignal;

global RawSignalRx1 RawSignalRx2 RawSignalRx3 RawSignalRx4 xTotalStationSync;


%scatter3(xTotalStationSync(1,:), xTotalStationSync(2,:), xTotalStationSync(3,:), 'filled');
%plot(RawSignalRx1(1,:)-RawSignalRx1(2,:))

% Tableau de 20x81920 -> 20 points avec 81920 samples pour chaque et
% une fréquence d'échantillonage de FsReference


pointsNum = length(RawSignalRx1(:,1));
xPos = zeros(1,pointsNum);
yPos = zeros(1,pointsNum);
tdaua = zeros(6,pointsNum);
tdoref = zeros(6,pointsNum);
tdoaStatTotale = zeros(6,pointsNum);


for i = 1:pointsNum
    [xPos(i),yPos(i),tdaua(:,i),tdoref(:,i)] = findpos(i);
    tdoaStatTotale(:,i) = trueTDOAGeom(xTotalStationSync(1,i),xTotalStationSync(2,i));
end


figure();
for i = 1:6
    plot(tdaua(i,:)*3*10^8);
    hold on;
end
figure();
for i = 1:6
    plot(tdoaStatTotale(i,:)*3*10^8);
    hold on;
end


% ---- PLOT POUR LE RAPPORT ----
%{
r1 = RawSignalRx1(1,:) - mean(RawSignalRx1(1,:));
r2 = RawSignalRx2(1,:) - mean(RawSignalRx2(1,:));
limiteTemp = 31000;
limiteLowRef = 36000;
limiteHighRef = 80000;
r1Balise = r1(1:limiteTemp);
r1Ref = r1(limiteLowRef:limiteHighRef);
r2Balise = r2(1:limiteTemp);
r2Ref = r2(limiteLowRef:limiteHighRef);
timeDelay = findDelay(r1Balise,r2Balise)
dist = timeDelay * 3 * 10^8
%}

% ---- PLOT ----
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
hold on;
scatter(xPos,yPos,'filled') % Position calcul�e
legend('Antennes','R�f�rence','Positions laser','Positions estim�es','Location','north')
xlabel('position [m]');
ylabel('position [m]');
grid on;


function [xPos,yPos,delayyy,tdoareff] = findpos(point)
    global RawSignalRx1 RawSignalRx2 RawSignalRx3 RawSignalRx4 xTotalStationSync;
    global xRef yRef tau;
    
    % On enlève la moyenne des signaux (partie DC)
    r1 = RawSignalRx1(point,:) - mean(RawSignalRx1(point,:));
    r2 = RawSignalRx2(point,:) - mean(RawSignalRx2(point,:));
    r3 = RawSignalRx3(point,:) - mean(RawSignalRx3(point,:));
    r4 = RawSignalRx4(point,:) - mean(RawSignalRx4(point,:));

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

    % Les TDOA qu'on devrait trouver
    TDOARefExact = trueTDOAGeom(xRef,yRef);

    % Delais des balises
    TDOABal = [findDelay(r1Balise,r2Balise),...
               findDelay(r1Balise,r3Balise),...
               findDelay(r1Balise,r4Balise),...
               findDelay(r2Balise,r3Balise),...
               findDelay(r2Balise,r4Balise),...
               findDelay(r3Balise,r4Balise)];

    % Delais de référence
    TDOARef = [findDelay(r1Ref,r2Ref),...
               findDelay(r1Ref,r3Ref),...
               findDelay(r1Ref,r4Ref),...
               findDelay(r2Ref,r3Ref),...
               findDelay(r2Ref,r4Ref),...
               findDelay(r3Ref,r4Ref)];

    % Correction des erreurs statiques
    errCorrection = TDOARef - TDOARefExact;

    % Delais corrigés des balises
    tau = TDOABal - errCorrection;
    
    x0 = [-3,1];
    [temp,~] = lsqnonlin(@func,x0,[],[],optimset('display','off'));
    xPos = temp(1);
    yPos = temp(2);
    delayyy = tau;
    tdoareff = errCorrection;
end

function [res] = upconvert(r)   
    R = upsample(r,3);
    L = length(R);
    spectre = fftshift(fft(R));    
    spectreFiltre = spectre;
    spectreFiltre(round(L/6):round(L*5/6)) = 0;
    res=ifft(ifftshift(spectreFiltre));
    
    %{
    global FsRawSignal;
    figure()
    Lr = length(r);
    fBB = (-Lr/2:Lr/2-1)*(FsRawSignal*10^-9/Lr);
    plot(fBB, abs(fftshift(fft(r))));
    title('Spectre en bande de base');
    xlabel('Fréquence [GHz]');
    figure()
    subplot(2,1,1);
    fRF = (-L/2:L/2-1)*(FsRawSignal*3*10^-9/L);
    plot(fRF, abs(spectre));
    title('Spectre en radiofréquences');
    xlabel('Fréquence [GHz]');
    subplot(2,1,2);
    plot(fRF, abs(spectreFiltre));
    title('Spectre en radiofréquences filtré');
    xlabel('Fréquence [GHz]');
    figure()
    plot(real(res))
    title('Signal temporel en radiofréquence');
    xlabel('Sample');
    %}
    
end

function timeDelay = findDelay(r1,r2)
    global Fs
    R1 = upconvert(r1);
    R2 = upconvert(r2);
    [acor,~] = xcorr(R1,R2);    
    % Maximum de cette corrélation -> index
    [~,maxIndex] = max(abs(acor));
    normax = maxIndex - length(acor)/2;
    timeDelay = normax/Fs;
    
    %{
    figure()
    L = length(acor);
    xAxe = -L/2:L/2-1;
    plot(xAxe,real(acor));
    title('Correlation');
    xlabel('Sample');
    grid on;
    %}
end


% Fonction qui retourne les TDOA théoriques sur base de la distance
% euclidienne entre les points
function TDOARef = trueTDOAGeom(x,y)
    global lightSpeed x1 x2 y1 y2;
    tau = zeros(1,6);
    tau(1) = (sqrt((x-x2(1))^2+(y-y2(1))^2)-sqrt((x-x1(1))^2+(y-y1(1))^2)) / lightSpeed;
    tau(2) = (sqrt((x-x2(2))^2+(y-y2(2))^2)-sqrt((x-x1(2))^2+(y-y1(2))^2)) / lightSpeed;
    tau(3) = (sqrt((x-x2(3))^2+(y-y2(3))^2)-sqrt((x-x1(3))^2+(y-y1(3))^2)) / lightSpeed;
    tau(4) = (sqrt((x-x2(4))^2+(y-y2(4))^2)-sqrt((x-x1(4))^2+(y-y1(4))^2)) / lightSpeed;
    tau(5) = (sqrt((x-x2(5))^2+(y-y2(5))^2)-sqrt((x-x1(5))^2+(y-y1(5))^2)) / lightSpeed;
    tau(6) = (sqrt((x-x2(6))^2+(y-y2(6))^2)-sqrt((x-x1(6))^2+(y-y1(6))^2)) / lightSpeed;
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

%{
figure()
scatter(1:length(tau),tau,'filled');
hold on;
scatter(1:length(TDOABal),TDOABal);
hold on;
scatter(1:length(errCorrection),errCorrection);
grid on;
legend('Corrected TDOA','Computed TDOA','Correction (TDOA in cables)');

figure()
scatter(1:length(errCorrection),errCorrection,'filled');
hold on;
scatter(1:length(TDOARef),TDOARef);
hold on;
scatter(1:length(TDOARefExact),TDOARefExact);
grid on;
legend('Correction','Computed TDOA','Geometrics TDOA');
%}

%{
% Fonction qui donne le TDOA entre les deux signaux.
% Calcule la corrélation puis le TDOA.
function timeDelay = findDelay(r1,r2)
    global Fs;
    
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
    timeDelay = sampleDelay/Fs;
    
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
%}
 
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


% Passe le signal vers de la radio fréquence.
% r  : signal temporel à mettre en radio-fréquences
% Fs : Fréquence du signal
% Fsample : Fréquence d'échantillonage
%{
function [r] = toRf(r, Fs)
    numZeros = 3;
    L_per = numZeros * length(r) - (numZeros-1);
    r_per = zeros(1,L_per);
    r_per(1:numZeros:end) = r;
    %R_per = abs(fftshift(fft(r_per))/L_per);
    f_per = numZeros*(-L_per/2:L_per/2-1)*(Fs/L_per);
    
    R_per = (fft(r_per));
    
    [~,indexFs] = min(abs(f_per-Fs));
    R_per(L_per-indexFs:indexFs) = 0;
    r = real(ifft(R_per));
    %f = f_per;
    %R = R_per;
end
%}
