function hfig = figure_size(SizeX,SizeY)
% Create a figure with size of SizeX and SizeY, in inch
% 5 inches by 3 inches by default
if nargin == 0
    SizeX=5;
    SizeY=3;
elseif nargin ~= 2
    error('invalid input argument');
end

hfig=figure('units','inches');
pos = get(hfig,'pos');
set(gcf,'pos',[pos(1) pos(2) SizeX SizeY]);

end