%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : X. Luo et al., Two-Dimensional Direction-of-Arrival Estimation Using Two Transform 
%               Matrices, IEEE, 2014.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 5 July 2017, Bandhit Suksiri,
%               Updated: 5 July 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out_ang_x_rad_vec, out_ang_z_rad_vec] = doa_frq_n_src ( ...
    src_mat, n_sen, d_sen, lamb, n_src, src_frq, cen_frq, ...
    ang_z_rad_vec)
    % n_sen_vec must be : [n_sen; 0; n_sen - 1]!
	is_eulr_ang = true;
    [n_tim, ~]  = size(src_mat);
    x_mat       = src_mat(:, 2: 1: n_sen + 1).';
    z_mat       = src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
    
    r_mat = zeros(n_sen, n_sen);
    for i = 1: 1: n_tim
        r_mat = r_mat + (z_mat(:, i) * x_mat(:, i)');
    end
    r_mat = r_mat / n_tim;
    
    [u_mat, d_mat, v_mat] = svd(r_mat);
    u_s_mat = u_mat(:, 1: n_src);
    d_s_mat = d_mat(1: n_src, 1: n_src);
    v_s_mat = v_mat(:, 1: n_src);
    u_n_mat = u_mat(:, n_src + 1: end);
%     d_n_mat = d_mat(n_src + 1: end, n_src + 1: end);
%     v_n_mat = v_mat(:, n_src + 1: end);
    
    n_ang_z     = size(ang_z_rad_vec, 1);
    norm_vec    = zeros(n_ang_z, 1);
    sqr_u_n_mat = u_n_mat * u_n_mat';
    for i = 1: 1: n_ang_z
        ang_rad_vec    = [ndef_ang, ndef_ang, ang_z_rad_vec(i, 1)];
        ster_vec       = l_fwd_ster_frq_fcn([0; 0; n_sen - 1], ...
                                            d_sen, ...
                                            ang_rad_vec, ...
                                            src_frq, ...
                                            cen_frq, ...
                                            lamb, ...
                                            is_eulr_ang);
        norm_vec(i, 1)  = 1 / (ster_vec' * sqr_u_n_mat * ster_vec);
    end
    norm_vec = abs(norm_vec);
    norm_vec = 10 .* log10(norm_vec);
    norm_vec(norm_vec <= 0) = 0;
    
    [pks_vec, locs_vec]   = findpeaks(norm_vec);
    [~, sort_pks_idx_vec] = sort(pks_vec, 'descend');
    n_pks                 = size(sort_pks_idx_vec, 1);
    if (n_pks < n_src)
        n_src = n_pks;
    end
    out_ang_z_rad_vec = ang_z_rad_vec(locs_vec(sort_pks_idx_vec(1 : n_src, 1), 1), 1);
    
    ang_rad_vec       = ndef_ang_mat(n_src, 3);
    ang_rad_vec(:, 3) = out_ang_z_rad_vec;
    a_z_mat           = l_fwd_ster_frq_fcn([0; 0; n_sen - 1], ...
                                           d_sen, ...
                                           ang_rad_vec, ...
                                           src_frq, ...
                                           cen_frq, ...
                                           lamb, ...
                                           is_eulr_ang);
    
	v_k_mat           = v_s_mat;
    d_k_mat           = d_s_mat;
    a_x_r_s_mat       = v_k_mat * d_k_mat * (pinv(a_z_mat) * u_s_mat)';
    a_1_x_r_s_mat     = a_x_r_s_mat(1: n_sen - 1, :);
    a_2_x_r_s_mat     = a_x_r_s_mat(2: n_sen,     :);
    
    out_ang_x_rad_vec = zeros(n_src, 1);
    for i = 1: 1: n_src
        out_ang_x_rad_vec(i, 1) = acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                                       angle(pinv(a_1_x_r_s_mat(:, i)) * a_2_x_r_s_mat(:, i)));
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%