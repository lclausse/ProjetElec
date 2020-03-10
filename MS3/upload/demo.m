clear all;
% assuming phi in [0, pi), theta in [0, 2*pi]
[PHI, THT] = meshgrid(linspace(0, pi*89/90, 90), linspace(0, 2*pi, 181));
gain_dB = 10*log10(abs(sin(THT)));
[PHI, THT, gain_dB] = patch_slot(PHI, THT, gain_dB);
[PHI, THT, gain_dB] = adjust_pattern_mesh_size(PHI, THT, gain_dB);
figure_size();
plot_antenna_pattern(gain_dB, PHI,THT, 20);
axis equal;


