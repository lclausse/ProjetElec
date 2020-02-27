clear
clc
importfile('Data_Measured.mat')

c = 299792458;
TDOA1 = -1.776939857891486e-08;
x1 = xReceivers(1,1)
y1 = xReceivers(2,1)
x2 = xReceivers(1,2)
y2 = xReceivers(2,2)

t = 5;

x = 2;
 
yArray = (t*((- t^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)*(- t^2 + 4*x^2 - 4*x*x1 - 4*x*x2 + x1^2 + 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))^(1/2) + t^2*y1 + t^2*y2 - x1^2*y1 + x1^2*y2 + x2^2*y1 - x2^2*y2 + y1*y2^2 + y1^2*y2 - y1^3 - y2^3 + 2*x*x1*y1 - 2*x*x1*y2 - 2*x*x2*y1 + 2*x*x2*y2)/(2*t^2 - 2*y1^2 + 4*y1*y2 - 2*y2^2)*(t^2*y1 - t*((- t^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)*(- t^2 + 4*x^2 - 4*x*x1 - 4*x*x2 + x1^2 + 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))^(1/2) + t^2*y2 - x1^2*y1 + x1^2*y2 + x2^2*y1 - x2^2*y2 + y1*y2^2 + y1^2*y2 - y1^3 - y2^3 + 2*x*x1*y1 - 2*x*x1*y2 - 2*x*x2*y1 + 2*x*x2*y2)/(2*t^2 - 2*y1^2 + 4*y1*y2 - 2*y2^2)


function importfile(fileToRead1)
    
    newData1 = load('-mat', fileToRead1);
    %Create new variables in the base workspace from those fields.
    vars = fieldnames(newData1);
    for i = 1:length(vars)
        assignin('base', vars{i}, newData1.(vars{i}));
    end
end