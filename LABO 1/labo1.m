clear all;
close all;
clc;
%importfile('Data_Measured.mat')
%importfile('Data_Synthetic.mat')
load('Data_Lab1_4.mat')

% Pour enlever le message d'erreur.
MSGID = 'MATLAB:declareGlobalBeforeUse';
warning('off', MSGID);
global lightSpeed;
lightSpeed = 299792458;

global RawSignalRx1 RawSignalRx2 RawSignalRx3 RawSignalRx4 xTotalStationSync;
global Fs Fsample FsRawSignal tau;

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


% Frequence d'echantillonnage
Fsample = FsReference;
% Frequence du signal
Fs = 3*FsRawSignal;

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

figure()
plot(tdaua(3,:));
hold on;
plot(tdoaStatTotale(3,:));
xlabel('Points');
ylabel('TDOA [s]');
legend('TDOA calcul�s','TDOA g�om�triques');
title('TDOA en radio-fr�quence');

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
% ---- END PLOT POUR LE RAPPORT ----


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
    global RawSignalRx1 RawSignalRx2 RawSignalRx3 RawSignalRx4;
    global xRef yRef tau;
    
    % On enlève la moyenne des signaux (partie DC)
    r1 = RawSignalRx1(point,:) - mean(RawSignalRx1(point,:));
    r2 = RawSignalRx2(point,:) - mean(RawSignalRx2(point,:));
    r3 = RawSignalRx3(point,:) - mean(RawSignalRx3(point,:));
    r4 = RawSignalRx4(point,:) - mean(RawSignalRx4(point,:));

    % Limites ou on coupe les signaux pour isoler partie de la balise et ref
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
    
    r1_Balise_RF = upconvert(r1Balise);
    r1_Ref_RF = upconvert(r1Ref);
    r2_Balise_RF = upconvert(r2Balise);
    r2_Ref_RF = upconvert(r2Ref);
    r3_Balise_RF = upconvert(r3Balise);
    r3_Ref_RF = upconvert(r3Ref);
    r4_Balise_RF = upconvert(r4Balise);
    r4_Ref_RF = upconvert(r4Ref);

    % Les TDOA qu'on devrait trouver
    TDOARefExact = trueTDOAGeom(xRef,yRef);

    % Delais des balises
    TDOABal = [findDelay(r1_Balise_RF,r2_Balise_RF),...
               findDelay(r1_Balise_RF,r3_Balise_RF),...
               findDelay(r1_Balise_RF,r4_Balise_RF),...
               findDelay(r2_Balise_RF,r3_Balise_RF),...
               findDelay(r2_Balise_RF,r4_Balise_RF),...
               findDelay(r3_Balise_RF,r4_Balise_RF)];

    % Delais de reference
    TDOARef = [findDelay(r1_Ref_RF,r2_Ref_RF),...
               findDelay(r1_Ref_RF,r3_Ref_RF),...
               findDelay(r1_Ref_RF,r4_Ref_RF),...
               findDelay(r2_Ref_RF,r3_Ref_RF),...
               findDelay(r2_Ref_RF,r4_Ref_RF),...
               findDelay(r3_Ref_RF,r4_Ref_RF)];
    
    %{
    % En bande de base
    TDOABal = [findDelay(r1Balise,r2Balise),...
               findDelay(r1Balise,r3Balise),...
               findDelay(r1Balise,r4Balise),...
               findDelay(r2Balise,r3Balise),...
               findDelay(r2Balise,r4Balise),...
               findDelay(r3Balise,r4Balise)];

    % Delais de r�f�rence
    TDOARef = [findDelay(r1Ref,r2Ref),...
               findDelay(r1Ref,r3Ref),...
               findDelay(r1Ref,r4Ref),...
               findDelay(r2Ref,r3Ref),...
               findDelay(r2Ref,r4Ref),...
               findDelay(r3Ref,r4Ref)];
    %}
    

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
    res = ifft(ifftshift(spectreFiltre));    
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


% IL FAUT LUI DONNER DES SIGNAUX EN RF MTN !!!
function timeDelay = findDelay(r1_RF,r2_RF)
    global Fs
    %r1_RF = upconvert(r1);
    %r2_RF = upconvert(r2);
    %[acor,~] = xcorr(r1_RF,r2_RF);    
    
    % Correlation fait maison
    e = fft(r1_RF);
    f = fft(r2_RF);
    acor = ifftshift(ifft(e.*conj(f)));
    
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


% Fonction qui retourne les TDOA theoriques sur base de la distance
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


% Fonction non lineaire a resoudre pour trouver les positions
function F = func(x)
    global x1 x2 y1 y2 lightSpeed tau;
    f1 = sqrt((x(1)-x2(1))^2+(x(2)-y2(1))^2)-sqrt((x(1)-x1(1))^2+(x(2)-y1(1))^2)-lightSpeed*tau(1); %A2-A1
    f2 = sqrt((x(1)-x2(2))^2+(x(2)-y2(2))^2)-sqrt((x(1)-x1(2))^2+(x(2)-y1(2))^2)-lightSpeed*tau(2); %A3-A1
    f3 = sqrt((x(1)-x2(3))^2+(x(2)-y2(3))^2)-sqrt((x(1)-x1(3))^2+(x(2)-y1(3))^2)-lightSpeed*tau(3); %A4-A1
    f4 = sqrt((x(1)-x2(4))^2+(x(2)-y2(4))^2)-sqrt((x(1)-x1(4))^2+(x(2)-y1(4))^2)-lightSpeed*tau(4); %A3-A2
    f5 = sqrt((x(1)-x2(5))^2+(x(2)-y2(5))^2)-sqrt((x(1)-x1(5))^2+(x(2)-y1(5))^2)-lightSpeed*tau(5); %A4-A2
    f6 = sqrt((x(1)-x2(6))^2+(x(2)-y2(6))^2)-sqrt((x(1)-x1(6))^2+(x(2)-y1(6))^2)-lightSpeed*tau(6); %A4-A3
    F = [f1, f2, f3, f4, f5, f6];
end






% ------------------------------------------------------------------------
% --------------------------- bullshit -----------------------------------
% ------------------------------------------------------------------------

 
% Passe le signal vers de la radio frequence.
% r  : signal temporel a mettre en radio-frequences
% Fs : Frequence du signal
% Fsample : Frequence d'echantillonage
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
