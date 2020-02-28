clear
clc
importfile('Data_Measured.mat')

global c;
c = 299792458;

testedValue = 90;

global x1;
x1 = [xReceivers(1,1), xReceivers(1,1), xReceivers(1,1), xReceivers(1,2), xReceivers(1,2),xReceivers(1,3)];
global y1;
y1 = [xReceivers(2,1), xReceivers(2,1), xReceivers(2,1), xReceivers(2,2), xReceivers(2,2),xReceivers(2,3)];
global x2;
x2 = [xReceivers(1,2), xReceivers(1,3), xReceivers(1,4), xReceivers(1,3), xReceivers(1,4),xReceivers(1,4)];
global y2;
y2 = [xReceivers(2,2), xReceivers(2,3), xReceivers(2,4), xReceivers(2,3), xReceivers(2,4),xReceivers(2,4)];


global tau;
tau = TDOA(:,testedValue);
x0 = [5,5];
xlsqn = lsqnonlin(@func,x0);
xEst = xlsqn(1);
yEst = xlsqn(2);


% ---- PLOT ----
scatter(xReceivers(1,:),xReceivers(2,:),'filled')
str = [" R1"," R2"," R3"," R4"];
text(xReceivers(1,:),xReceivers(2,:), str)
hold on;

scatter(xEst,yEst,'*')
hold on;

scatter(xTag(1,testedValue),xTag(2,testedValue),'o')

syms xSym ySym
x = [xSym, ySym];
f1 =  @(xSym, ySym) sqrt((xSym-x2(1)).^2+(ySym-y2(1)).^2)-sqrt((xSym-x1(1)).^2+(ySym-y1(1)).^2) - c*tau(1);
f2 =  @(xSym, ySym) sqrt((xSym-x2(2)).^2+(ySym-y2(2)).^2)-sqrt((xSym-x1(2)).^2+(ySym-y1(2)).^2) - c*tau(2);
f3 =  @(xSym, ySym) sqrt((xSym-x2(3)).^2+(ySym-y2(3)).^2)-sqrt((xSym-x1(3)).^2+(ySym-y1(3)).^2) - c*tau(3);
f4 =  @(xSym, ySym) sqrt((xSym-x2(4)).^2+(ySym-y2(4)).^2)-sqrt((xSym-x1(4)).^2+(ySym-y1(4)).^2) - c*tau(4);
f5 =  @(xSym, ySym) sqrt((xSym-x2(5)).^2+(ySym-y2(5)).^2)-sqrt((xSym-x1(5)).^2+(ySym-y1(5)).^2) - c*tau(5);
f6 =  @(xSym, ySym) sqrt((xSym-x2(6)).^2+(ySym-y2(6)).^2)-sqrt((xSym-x1(6)).^2+(ySym-y1(6)).^2) - c*tau(6);

fimplicit(f1,[0 10 -1 9])
fimplicit(f2,[0 10 -1 9])
fimplicit(f3,[0 10 -1 9])
fimplicit(f4,[0 10 -1 9])
fimplicit(f5,[0 10 -1 9])
fimplicit(f6,[0 10 -1 9])

legend('Récepteurs','Estimation','True pos')


function F = func(x)
    global x1 x2 y1 y2 c tau;
    xSym = x(1);
    ySym = x(2);
    
    f1 = sqrt((xSym-x2(1)).^2+(ySym-y2(1)).^2)-sqrt((xSym-x1(1)).^2+(ySym-y1(1)).^2) - c*tau(1);
    f2 = sqrt((xSym-x2(2)).^2+(ySym-y2(2)).^2)-sqrt((xSym-x1(2)).^2+(ySym-y1(2)).^2) - c*tau(2);
    f3 = sqrt((xSym-x2(3)).^2+(ySym-y2(3)).^2)-sqrt((xSym-x1(3)).^2+(ySym-y1(3)).^2) - c*tau(3);
    f4 = sqrt((xSym-x2(4)).^2+(ySym-y2(4)).^2)-sqrt((xSym-x1(4)).^2+(ySym-y1(4)).^2) - c*tau(4);
    f5 = sqrt((xSym-x2(5)).^2+(ySym-y2(5)).^2)-sqrt((xSym-x1(5)).^2+(ySym-y1(5)).^2) - c*tau(5);
    f6 = sqrt((xSym-x2(6)).^2+(ySym-y2(6)).^2)-sqrt((xSym-x1(6)).^2+(ySym-y1(6)).^2) - c*tau(6);
    
    F = [f1,f2,f3,f4,f5,f6];
end



function importfile(fileToRead1)
    newData1 = load('-mat', fileToRead1);
    %Create new variables in the base workspace from those fields.
    vars = fieldnames(newData1);
    for i = 1:length(vars)
        assignin('base', vars{i}, newData1.(vars{i}));
    end
end