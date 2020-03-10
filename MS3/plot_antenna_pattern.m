function handle=plot_antenna_pattern(gain,PHI,THETA,displacement)
% 'displacement' determines the smallest gain shown in the plot, anything
% smaller than this value will be forcely evaluated to this value

color_scale=10;
color_max=color_scale*(ceil(max(max(gain)))+displacement);
N_ticks=5;
    

% coordinate transform
gain=color_scale*(gain+displacement);
gain(gain<0)=0;
X=gain.*sin(THETA).*cos(PHI);
Y=gain.*sin(THETA).*sin(PHI);
Z=gain.*cos(THETA);
% plot
handle=warp(X,Y,Z,gain,jet(color_max));
% set up a colorbar
ticks=linspace(0,color_max,N_ticks);
ticklabels=textscan(num2str(round(ticks/color_scale-displacement)),'%s');
ticklabels=ticklabels{1,1}';
ticklabels{1,N_ticks}=strcat(ticklabels{1,N_ticks}); % ,'(dBi)'
colorbar('Ticks',ticks,'TickLabels',ticklabels);

set(gca,'Visible','off');
shading faceted;

end

