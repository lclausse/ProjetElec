clear
clc
importfile('Data_Measured.mat')
%importfile('Data_Synthetic.mat')
%importfile('Data_Lab1_1.mat')


global lightSpeed;
lightSpeed = 299792458;

% xTag : Positions au laser
% xReceivers : Positions des balises au laser

global x1;
x1 = [xReceivers(1,1), xReceivers(1,1), xReceivers(1,1), xReceivers(1,2), xReceivers(1,2),xReceivers(1,3)];
global y1;
y1 = [xReceivers(2,1), xReceivers(2,1), xReceivers(2,1), xReceivers(2,2), xReceivers(2,2),xReceivers(2,3)];
global x2;
x2 = [xReceivers(1,2), xReceivers(1,3), xReceivers(1,4), xReceivers(1,3), xReceivers(1,4),xReceivers(1,4)];
global y2;
y2 = [xReceivers(2,2), xReceivers(2,3), xReceivers(2,4), xReceivers(2,3), xReceivers(2,4),xReceivers(2,4)];

xEst = zeros(1,length(TDOA));
yEst = zeros(1,length(TDOA));
resnorm = zeros(1,length(TDOA));
global tau;


for i = 1:length(TDOA)
    tau = TDOA(:,i);
    x0 = [5,5];
    [temp,resnorm(i)] = lsqnonlin(@func,x0);
    xEst(i) = temp(1);
    yEst(i) = temp(2);
end


l = xTag(1,:);
err = immse(xEst,l);
% ---- PLOT ----
scatter(xReceivers(1,:),xReceivers(2,:),'filled')
str = [" R1"," R2"," R3"," R4"];
text(xReceivers(1,:),xReceivers(2,:), str)
hold on;
scatter(xEst,yEst,'*')
hold on;
text(-8, 12, "MSE = " + err)
hold on;
scatter(xTag(1,:),xTag(2,:),'k','.')
legend('Récepteurs','Estimation','True pos')


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



function importfile(fileToRead1)
    newData1 = load('-mat', fileToRead1);
    %Create new variables in the base workspace from those fields.
    vars = fieldnames(newData1);
    for i = 1:length(vars)
        assignin('base', vars{i}, newData1.(vars{i}));
    end
end