%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Abstract Forward Equation of Steering Vector
%
% Description : Abstract Forward Equation of Steering Vector
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 4 June 2017, Bandhit Suksiri,
%               Updated: 4 June 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ster_mat = abs_fwd_ster_fcn (pos_mat, ang_rad_mat, src_frq, c, is_eulr_ang)
    if is_eulr_ang == true
        ang_x_rad_vec = ang_rad_mat(:, 1);
        ang_y_rad_vec = ang_rad_mat(:, 2);
        ang_z_rad_vec = ang_rad_mat(:, 3);
        ster_mat      = ((2 * pi * src_frq * 1i) / c) .* ...
                        ((pos_mat(:, 1) * cos(ang_x_rad_vec')) + ...
                         (pos_mat(:, 2) * cos(ang_y_rad_vec')) + ...
                         (pos_mat(:, 3) * cos(ang_z_rad_vec')));
    else
        az_rad_vec = ang_rad_mat(:, 1);
        el_rad_vec = ang_rad_mat(:, 2);
        ster_mat   = ((2 * pi * src_frq * 1i) / c) .* ...
                     ((pos_mat(:, 1) * (sin(el_rad_vec) .* cos(az_rad_vec))') + ...
                      (pos_mat(:, 2) * (sin(el_rad_vec) .* sin(az_rad_vec))') + ...
                      (pos_mat(:, 3) * cos(el_rad_vec)'));
    end
    ster_mat = exp(ster_mat);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%