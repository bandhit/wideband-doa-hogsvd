%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : J. F. Gu et al., Joint SVD of Two Cross-Correlation Matrices to Achieve Automatic 
%               Pairing in 2-D Angle Estimation Problems, IEEE, 2007.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 7 July 2017, Bandhit Suksiri,
%               Updated: 7 July 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out_ang_x_rad_vec, out_ang_z_rad_vec] = doa_frq_n_src ( ...
    src_mat, n_sen, d_sen, lamb, n_src, src_frq, cen_frq, ...
    ang_x_rad_vec)
    % n_sen_vec must be : [n_sen; 0; n_sen - 1]!
	is_eulr_ang = true;
    [n_tim, ~]  = size(src_mat);
    x_mat       = src_mat(:, 2: 1: n_sen + 1).';
    z_mat       = src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
    
    r_mat = zeros(n_sen, n_sen);
    for i = 1: 1: n_tim
        r_mat = r_mat + (x_mat(:, i) * z_mat(:, i)');
    end
    r_mat = r_mat / n_tim;
    
    r_1_mat = r_mat(:, 1: n_sen - 1);
    r_2_mat = r_mat(:, 2: n_sen);
    new_r_mat = [r_1_mat; r_2_mat];
    
    [u_mat, ~, ~] = svd(new_r_mat);
    u_1_mat       = u_mat  (                   :, 1: n_src);
    u_1_1_mat     = u_1_mat(        1:     n_sen, :);
    u_1_2_mat     = u_1_mat(n_sen + 1: 2 * n_sen, :);
    
    [inv_t_mat, ohm_mat] = eig(pinv(u_1_1_mat) * u_1_2_mat);
    
    a_x_mat           = u_1_1_mat * inv_t_mat;
    n_ang_x           = size(ang_x_rad_vec, 1);
    out_ang_x_rad_vec = zeros(n_src, 1);
    for i = 1: 1: n_src
        obj_vec = zeros(n_ang_x, 1);
        for j = 1: 1: n_ang_x
            ang_rad_vec   = [ang_x_rad_vec(j, 1), ndef_ang, ndef_ang];
            ster_vec      = l_fwd_ster_frq_fcn([n_sen; 0; 0], ...
                                               d_sen, ...
                                               ang_rad_vec, ...
                                               src_frq, ...
                                               cen_frq, ...
                                               lamb, ...
                                               is_eulr_ang);
            ster_vec      = ster_vec(2: end, :);
            obj_vec(j, 1) = abs(ster_vec' * a_x_mat(:, i)) / (norm(ster_vec, 2) ^ 2);
        end
        [~, idx] = max(obj_vec);
        out_ang_x_rad_vec(i, 1) = ang_x_rad_vec(idx, 1);
    end
    
    out_ang_z_rad_vec = acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                             angle(diag(ohm_mat)));
    out_ang_z_rad_vec = pi - out_ang_z_rad_vec;
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%