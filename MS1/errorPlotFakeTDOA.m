clear
clc
close all

load('Data_Measured.mat')

global c tau;
c = 299792458;

testedValue = 40;

global x1;
x1 = [xReceivers(1,1), xReceivers(1,1), xReceivers(1,1), xReceivers(1,2), xReceivers(1,2),xReceivers(1,3)];
global y1;
y1 = [xReceivers(2,1), xReceivers(2,1), xReceivers(2,1), xReceivers(2,2), xReceivers(2,2),xReceivers(2,3)];
global x2;
x2 = [xReceivers(1,2), xReceivers(1,3), xReceivers(1,4), xReceivers(1,3), xReceivers(1,4),xReceivers(1,4)];
global y2;
y2 = [xReceivers(2,2), xReceivers(2,3), xReceivers(2,4), xReceivers(2,3), xReceivers(2,4),xReceivers(2,4)];

% Méthode :
% On crée un maillage géométrique
% On calcule les TDOA exacts
% On met une erreur sur ces TDOA
% On calcule la position trouvée avec ces faux tdoa
% On prend l'erreur de lsqnonlin
% On plot sur le mesh l'erreur

TDOARef = trueTDOARef(1,1)

x = -10:1:18;
y = -10:1:18;
err = zeros(length(x),length(y));

opts = optimset('Display','off');

for i = 1:length(x)
    for j = 1:length(y)
        TDOARef = trueTDOARef(x(i),y(j));
        error = 2 * 10^-10; % rand(1,6) * 0.5 * 10^-9;
        tau = TDOARef + error;

        x0 = [5,5];
        [temp,err(i,j)] = lsqnonlin(@func, x0, [], [], opts);
        err(i,j) = sqrt((x(i)-temp(1))^2 + (y(j)-temp(2))^2);
    end
end

% Retirer gros piques
%[indX,indY] = find(err > 6);
%err(indX,indY) = 6;

% ---- PLOT Error 3D ----

figure()
surf(x,y,err);
hold on;
height = 1.1;% * max(max(err));
scatter3(xReceivers(1,:),xReceivers(2,:),zeros(1,length(xReceivers(2,:)))+height,'filled','r')
hold on;
for i = 1:length(xReceivers(1,:))
    plot3([xReceivers(1,i) xReceivers(1,i)], [xReceivers(2,i) xReceivers(2,i)], [zeros(1,length(xReceivers(2,i)))+height zeros(1,length(xReceivers(2,i)))],'r');
    hold on;
end
xlabel('position [m]');
ylabel('position [m]');
zlabel('Erreur [m]');
title('Erreur de position');

% ---- END PLOT Error 3D ----


% ---- PLOT Error 2D ----
figure()
pcolor(x,y,log10(err));
shading interp;
hold on;
scatter(xReceivers(1,:),xReceivers(2,:),'filled','r')
xlabel('position [m]');
ylabel('position [m]');
title('Erreur de position');
legend('Error','Récepteurs');
% ---- END PLOT Error 2D ----


% ---- PLOT ----
%{
figure();
scatter(xReceivers(1,:),xReceivers(2,:),'filled')
str = [" R1"," R2"," R3"," R4"];
text(xReceivers(1,:),xReceivers(2,:), str)
hold on;
scatter(xEst,yEst,'k*')
hold on;
scatter(xTag(1,testedValue),xTag(2,testedValue),'ko')

syms xSym ySym
x = [xSym, ySym];
f1 =  @(xSym, ySym) sqrt((xSym-x2(1)).^2+(ySym-y2(1)).^2)-sqrt((xSym-x1(1)).^2+(ySym-y1(1)).^2) - c*tau(1);
f2 =  @(xSym, ySym) sqrt((xSym-x2(2)).^2+(ySym-y2(2)).^2)-sqrt((xSym-x1(2)).^2+(ySym-y1(2)).^2) - c*tau(2);
f3 =  @(xSym, ySym) sqrt((xSym-x2(3)).^2+(ySym-y2(3)).^2)-sqrt((xSym-x1(3)).^2+(ySym-y1(3)).^2) - c*tau(3);
f4 =  @(xSym, ySym) sqrt((xSym-x2(4)).^2+(ySym-y2(4)).^2)-sqrt((xSym-x1(4)).^2+(ySym-y1(4)).^2) - c*tau(4);
f5 =  @(xSym, ySym) sqrt((xSym-x2(5)).^2+(ySym-y2(5)).^2)-sqrt((xSym-x1(5)).^2+(ySym-y1(5)).^2) - c*tau(5);
f6 =  @(xSym, ySym) sqrt((xSym-x2(6)).^2+(ySym-y2(6)).^2)-sqrt((xSym-x1(6)).^2+(ySym-y1(6)).^2) - c*tau(6);

limits = [0 10 -1 9];
fimplicit(f1,limits,'r')
fimplicit(f2,limits,'r')
fimplicit(f3,limits,'r')
fimplicit(f4,limits,'r')
fimplicit(f5,limits,'r')
fimplicit(f6,limits,'r')

legend('Récepteurs','Estimation','True pos')
grid on;
%}
% ---- END PLOT ----






function F = func(x)
    global x1 x2 y1 y2 c tau;
    f1 = sqrt((x(1)-x2(1))^2+(x(2)-y2(1))^2)-sqrt((x(1)-x1(1))^2+(x(2)-y1(1))^2)-c*tau(1);
    f2 = sqrt((x(1)-x2(2))^2+(x(2)-y2(2))^2)-sqrt((x(1)-x1(2))^2+(x(2)-y1(2))^2)-c*tau(2);
    f3 = sqrt((x(1)-x2(3))^2+(x(2)-y2(3))^2)-sqrt((x(1)-x1(3))^2+(x(2)-y1(3))^2)-c*tau(3);
    f4 = sqrt((x(1)-x2(4))^2+(x(2)-y2(4))^2)-sqrt((x(1)-x1(4))^2+(x(2)-y1(4))^2)-c*tau(4);
    f5 = sqrt((x(1)-x2(5))^2+(x(2)-y2(5))^2)-sqrt((x(1)-x1(5))^2+(x(2)-y1(5))^2)-c*tau(5);
    f6 = sqrt((x(1)-x2(6))^2+(x(2)-y2(6))^2)-sqrt((x(1)-x1(6))^2+(x(2)-y1(6))^2)-c*tau(6);
    F = [f1, f2, f3, f4, f5, f6];
end

function TDOARef = trueTDOARef(x,y)
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