clear
clc

importfile('Data.mat')


function importfile(fileToRead1)
    %  Imports data from the specified file
    newData1 = load('-mat', fileToRead1);

    % Create new variables in the base workspace from those fields.
    vars = fieldnames(newData1);
    for i = 1:length(vars)
        assignin('base', vars{i}, newData1.(vars{i}));
    end
end