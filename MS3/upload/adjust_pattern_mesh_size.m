function [ PHI_new,THT_new,gain_new ] = adjust_pattern_mesh_size(PHI,THT,gain)
%ADJUST_PATTERN_MESH_SIZE Summary of this function goes here
%   Detailed explanation goes here
N=30;
PHI_max=max(PHI(1,:));
PHI_min=min(PHI(1,:));
THT_max=max(THT(:,1));
THT_min=min(THT(:,1));
R_phi_tht=(PHI_max-PHI_min)/(THT_max-THT_min); % ratio of PHI range to THT range
if R_phi_tht > 1
    Nphi=round(N*R_phi_tht);
    Ntht=N;
else
    Ntht=round(N/R_phi_tht);
    Nphi=N;
end
[PHI_new,THT_new]=meshgrid(linspace(PHI_min,PHI_max,Nphi),...
    linspace(THT_min,THT_max,Ntht));
gain_new=interp2(PHI,THT,gain,PHI_new,THT_new);

end

