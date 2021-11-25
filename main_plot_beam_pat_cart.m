% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;
%clc;

% include
addpath('miscellaneous/');
addpath('model/');
addpath('third party/');

% configuration
n_mic          = 4;
n_mic_vec      = [n_mic; 0; n_mic - 1];
c              = 340;
cen_frq        = 3.4e3;
lamb           = c / cen_frq;
d_mic          = lamb / 2;
az_ang_rad_vec = (-pi: 1e-2: pi)';
el_ang_rad_vec = (-pi: 1e-2: pi)';
az_ang_deg_vec = az_ang_rad_vec * 180 / pi;
el_ang_deg_vec = el_ang_rad_vec * 180 / pi;

figure;
subplot(2, 1, 1);
beam_surf_plot(el_ang_deg_vec, az_ang_deg_vec, ...
    baem_pat_fcn(az_ang_rad_vec, el_ang_rad_vec, n_mic_vec, d_mic, cen_frq / 1, cen_frq, lamb));
subplot(2, 1, 2);
beam_surf_plot(el_ang_deg_vec, az_ang_deg_vec, ...
    baem_pat_fcn(az_ang_rad_vec, el_ang_rad_vec, n_mic_vec, d_mic, cen_frq / 4, cen_frq, lamb));

function af_mat = baem_pat_fcn (az_ang_rad_vec, el_ang_rad_vec, n_sen_vec, d_sen, src_frq, cen_frq, lamb)
    n_az   = size(az_ang_rad_vec, 1);
    n_el   = size(el_ang_rad_vec, 1);
    af_mat = zeros(n_az, n_el);
    for i = 1: 1: n_az
        for j = 1: 1: n_el
            af_mat(i, j) = mean(l_fwd_ster_frq_fcn(n_sen_vec, d_sen, ...
                                                  [az_ang_rad_vec(i, 1), el_ang_rad_vec(j, 1)], ...
                                                  src_frq, cen_frq, lamb, false));
        end
    end
    af_mat = abs(af_mat);
end

function beam_surf_plot (el_ang_deg_vec, az_ang_deg_vec, af_mat)
    surf(el_ang_deg_vec, az_ang_deg_vec, af_mat, ...
         'EdgeColor', 'none', 'LineStyle', 'none');
    xlabel('Elevation, Degree');
    ylabel('Azimuth, Degree');
    shading('interp');
    colormap(jet(1024));
    bar = colorbar;
    ylabel(bar, 'Magnitude');
    view(-45, +45);
    xlim([el_ang_deg_vec(1, 1), el_ang_deg_vec(end, 1)]);
    ylim([az_ang_deg_vec(1, 1), az_ang_deg_vec(end, 1)]);
    zlim([0, 1]);
    caxis([0, 1]);
    grid on;
    grid minor;
end
