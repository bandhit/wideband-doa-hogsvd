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
c              = 343;
cen_frq        = 2.0e3;
lamb           = c / cen_frq;
d_mic          = lamb / 2;
az_ang_rad_vec = (-pi: 1e-2: pi).';
el_ang_rad_vec = (-pi: 1e-2: pi).';

fig_obj   = figure;
file_name = 'beam_pat_sph';
x_val     = 1 * 21.0 / 1;
y_val     = 2 * 29.7 / 3;

subplot(2, 2, 1);
[af_mat, az_ang_rad_mat, el_ang_rad_mat] = baem_pat_fcn(az_ang_rad_vec, el_ang_rad_vec, n_mic_vec, d_mic, cen_frq / 2, cen_frq, lamb);
beam_surf_sph_plot(af_mat, az_ang_rad_mat, el_ang_rad_mat);
subplot(2, 2, 2);
[af_mat, az_ang_rad_mat, el_ang_rad_mat] = baem_pat_fcn(az_ang_rad_vec, el_ang_rad_vec, n_mic_vec, d_mic, cen_frq * 1, cen_frq, lamb);
beam_surf_sph_plot(af_mat, az_ang_rad_mat, el_ang_rad_mat);
subplot(2, 2, 3);
[af_mat, az_ang_rad_mat, el_ang_rad_mat] = baem_pat_fcn(az_ang_rad_vec, el_ang_rad_vec, n_mic_vec, d_mic, cen_frq * 2, cen_frq, lamb);
beam_surf_sph_plot(af_mat, az_ang_rad_mat, el_ang_rad_mat);
subplot(2, 2, 4);
[af_mat, az_ang_rad_mat, el_ang_rad_mat] = baem_pat_fcn(az_ang_rad_vec, el_ang_rad_vec, n_mic_vec, d_mic, cen_frq * 4, cen_frq, lamb);
beam_surf_sph_plot(af_mat, az_ang_rad_mat, el_ang_rad_mat);

% save_pdf(file_name, fig_obj, x_val, y_val);
% save_png(file_name, fig_obj, x_val, y_val, '600');
% pause(1);
% close(fig_obj);

function [af_mat, az_ang_rad_mat, el_ang_rad_mat] = baem_pat_fcn (az_ang_rad_vec, el_ang_rad_vec, n_sen_vec, d_sen, src_frq, cen_frq, lamb)
    n_az           = size(az_ang_rad_vec, 1);
    n_el           = size(el_ang_rad_vec, 1);
    af_mat         = zeros(n_el, n_az);
    az_ang_rad_mat = zeros(n_el, n_az);
    el_ang_rad_mat = zeros(n_el, n_az);
    for i = 1: 1: n_az
        for j = 1: 1: n_el
            af_mat(j, i) = mean(l_fwd_ster_frq_fcn(n_sen_vec, d_sen, ...
                                                   [az_ang_rad_vec(i, 1), el_ang_rad_vec(j, 1)], ...
                                                   src_frq, cen_frq, lamb, false));
            az_ang_rad_mat(j, i) = az_ang_rad_vec(i, 1);
            el_ang_rad_mat(j, i) = el_ang_rad_vec(j, 1);
        end
    end
    af_mat = abs(af_mat);
end

function beam_surf_sph_plot (af_mat, az_ang_rad_mat, el_ang_rad_mat)
    [x_mat, y_mat, z_mat] = sph_to_cart_rad(af_mat, az_ang_rad_mat, el_ang_rad_mat + (pi / 2));
    surf(x_mat, y_mat, z_mat, af_mat, ...
         'EdgeColor', 'none', 'LineStyle', 'none');
    xlabel('x, Coordinate');
    ylabel('y, Coordinate');
    zlabel('z, Coordinate');
    shading('interp');
    colormap(jet(1024));
    bar = colorbar('Location', 'northoutside');
    ylabel(bar, 'Magnitude');
    xlim([-1, +1]);
    ylim([-1, +1]);
    zlim([-1, +1]);
    view(-20, +20);
    caxis([0, 1]);
    grid on;
    grid minor;
end

function [x_mat, y_mat, z_mat] = sph_to_cart_deg(af_mat, az_ang_deg_mat, el_ang_deg_mat)
    x_mat = af_mat .* cosd(el_ang_deg_mat) .* cosd(az_ang_deg_mat);
    y_mat = af_mat .* cosd(el_ang_deg_mat) .* sind(az_ang_deg_mat);
    z_mat = af_mat .* sind(el_ang_deg_mat);
end

function [x_mat, y_mat, z_mat] = sph_to_cart_rad(af_mat, az_ang_rad_mat, el_ang_rad_mat)
    x_mat = af_mat .* cos(el_ang_rad_mat) .* cos(az_ang_rad_mat);
    y_mat = af_mat .* cos(el_ang_rad_mat) .* sin(az_ang_rad_mat);
    z_mat = af_mat .* sin(el_ang_rad_mat);
end

