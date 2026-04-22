% -----------------------------
% Define your plane
% -----------------------------
clear all;
close all;
%if you have a gsd file, you can use uncomment and use the functioin below
%to extract your frame(s).
% system=Extract_GSD('SPRING_ONLY_X_TO_Y_LxLyLz_Make_PCNDNano81_r2_N18NA9_nG110_VISUALIZATION.gsd','steps',-1);
system = Extract_XML('LxLyLz_Make_PCNDNano75_r2_N18NA9_nG110_RESIZED_BOXBEGIN_THERMO_INT.xml');

%Specifiy total length of polymer arm
arm_length=18;
%This specifies L assuming half beads are C and half are D
L=arm_length/2;
%Specify along which dimension the lamlellae are oriented
dimension_alignment = 'X';
D_ind=find(system.attype=='D');
C_ind=find(system.attype=='C');
D_interface=system.pos(D_ind(1:L:end-L),:);

if dimension_alignment == 'X'
    dim_col = 1;
    dim_size = system.dim(1);
elseif dimension_alignment == 'Y'
    dim_col = 2;
    dim_size = system.dim(2);
elseif dimension_alignment == 'Z'
    dim_col = 3;
    dim_size = system.dim(3);
end

num_bins = 100;
x_edges = linspace(-dim_size/2, dim_size/2, num_bins+1);
x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
[~, peak_locs_D] = findpeaks(histcounts(D_interface(:,dim_col), x_edges), x_centers, 'MinPeakProminence', max(histcounts(D_interface(:,dim_col), x_edges))*0.3);
% Midpoints between peaks are the natural dividers
boundaries = (peak_locs_D(1:end-1) + peak_locs_D(2:end)) / 2;
% Add -inf and +inf as outer boundaries
boundaries = [-inf, boundaries, inf];
% Assign each interface
D_lam_pos_1 = D_interface(D_interface(:,dim_col) > boundaries(1) & D_interface(:,dim_col) < boundaries(2), :);
D_lam_pos_2 = D_interface(D_interface(:,dim_col) > boundaries(2) & D_interface(:,dim_col) < boundaries(3), :);
D_lam_pos_3 = D_interface(D_interface(:,dim_col) > boundaries(3) & D_interface(:,dim_col) < boundaries(4), :);
D_lam_pos_4 = D_interface(D_interface(:,dim_col) > boundaries(4) & D_interface(:,dim_col) < boundaries(5), :);
D_lam_pos_5 = D_interface(D_interface(:,dim_col) > boundaries(5) & D_interface(:,dim_col) < boundaries(6), :);
D_lam_pos_6 = D_interface(D_interface(:,dim_col) > boundaries(6) & D_interface(:,dim_col) < boundaries(7), :);
total_captured = size(D_lam_pos_1,1) + size(D_lam_pos_2,1) + size(D_lam_pos_3,1) + ...
                 size(D_lam_pos_4,1) + size(D_lam_pos_5,1) + size(D_lam_pos_6,1);
fprintf('Total D_interface beads: %d\n', size(D_interface,1));
fprintf('Total captured: %d\n', total_captured);
fprintf('Missing: %d\n', size(D_interface,1) - total_captured);
D_lam_pos_all = {D_lam_pos_1, D_lam_pos_2, D_lam_pos_3, D_lam_pos_4, D_lam_pos_5, D_lam_pos_6};
A_total = 0;
A_all = zeros(1,6);
for surf = 1:6
    D = D_lam_pos_all{surf};
    x = D(:,1);
    y = D(:,2);
    z = D(:,3);
    if dimension_alignment == 'X'
        tri = delaunay(y, z);
        area_of_plane = system.dim(2)*system.dim(3);
    elseif dimension_alignment == 'Y'
        tri = delaunay(x,z);
        area_of_plane = system.dim(1)*system.dim(3);
    elseif dimension_alignment=='Z'
        tri = delaunay(x,y);
        area_of_plane = system.dim(1)*system.dim(2);
    else
        disp('Dimension not properly specified')
    end
    A = 0;
    for i = 1:size(tri,1)
        p1 = D(tri(i,1),:);
        p2 = D(tri(i,2),:);
        p3 = D(tri(i,3),:);
        A = A + 0.5 * norm(cross(p2 - p1, p3 - p1));
    end
    A_all(surf) = A;
    A_total = A_total + A;
    fprintf('Interface %d: %.4f\n', surf, A);
end
A_avg = A_total / 6;
A_per_V=sum(A_all)/(system.dim(1)*system.dim(2)*system.dim(3));
fprintf('Area_of_rough_surface,Area_of_smooth_surface,Ratio_Rough/Smooth,Interface_Area_per_Volume\n');
fprintf('%.4f,%.4f,%.4f,%.4f\n', A_avg, area_of_plane, A_avg/area_of_plane, A_per_V);
