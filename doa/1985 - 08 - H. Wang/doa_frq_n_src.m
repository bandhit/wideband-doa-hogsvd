%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : H. Wang et al., Coherent Signal-subspace Processing for the Detection and 
%               Estimation of Angles of Arrival of Multiple Wide-band Sources, IEEE, 1985.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 11 September 2017, Bandhit Suksiri,
%               Updated: 11 September 2017, Bandhit Suksiri.
%
% Copyright 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [init_ang_x_rad_vec, init_ang_z_rad_vec] = ...
    doa_frq_n_src (src_mat3, frq_vec, n_sen, d_sen, lamb, n_src, cen_frq, ...
    pre_ang_x_rad_vec, pre_ang_z_rad_vec, n_iter, ang_rad_vec)
    [n_frq, n_tim, ~]  = size(src_mat3);
    n_sen_vec          = [n_sen; 0; n_sen - 1];
    is_eulr_ang        = true;
    init_ang_x_rad_vec = pre_ang_x_rad_vec;
    init_ang_z_rad_vec = pre_ang_z_rad_vec;
    init_n_src_x       = size(init_ang_x_rad_vec, 1);
    init_n_src_z       = size(init_ang_z_rad_vec, 1);
    
    for i_iter = 1: 1: n_iter
        init_ang_x_rad_mat = [init_ang_x_rad_vec, zeros(init_n_src_x, 1), zeros(init_n_src_x, 1)];
        a_x_fc_mat         = l_fwd_ster_fcn(n_sen_vec, d_sen, init_ang_x_rad_mat, lamb, is_eulr_ang);
        a_x_fc_mat         = a_x_fc_mat(2: 1: n_sen + 1, :);
        init_ang_z_rad_mat = [zeros(init_n_src_z, 1), zeros(init_n_src_z, 1), init_ang_z_rad_vec];
        a_z_fc_mat         = l_fwd_ster_fcn(n_sen_vec, d_sen, init_ang_z_rad_mat, lamb, is_eulr_ang);
        a_z_fc_mat         = a_z_fc_mat([1, n_sen + 2: 1: 2 * n_sen], :);
        b_x_fc_mat         = [zeros(init_n_src_x, n_sen - init_n_src_x);
                              eye(n_sen - init_n_src_x)];
        b_x_fi_mat         = b_x_fc_mat;
        b_z_fc_mat         = [zeros(init_n_src_z, n_sen - init_n_src_z);
                              eye(n_sen - init_n_src_z)];
        b_z_fi_mat         = b_z_fc_mat;

        r_x_mat = zeros(n_sen, n_sen);
        r_z_mat = zeros(n_sen, n_sen);
        for i_frq = 1: 1: n_frq
            src_frq    = frq_vec(i_frq, 1);
            fi_src_mat = squeeze(src_mat3(i_frq, :, :));
            fi_x_mat   = fi_src_mat(:, 2: 1: n_sen + 1).';
            fi_z_mat   = fi_src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';

            fi_x_fi_x_r_mat = zeros(n_sen, n_sen);
            fi_z_fi_z_r_mat = zeros(n_sen, n_sen);
            for i_tim = 1: 1: n_tim
                fi_x_fi_x_r_mat = fi_x_fi_x_r_mat + (fi_x_mat(:, i_tim) * fi_x_mat(:, i_tim)');
                fi_z_fi_z_r_mat = fi_z_fi_z_r_mat + (fi_z_mat(:, i_tim) * fi_z_mat(:, i_tim)');
            end
            fi_x_fi_x_r_mat = fi_x_fi_x_r_mat / n_tim;
            fi_z_fi_z_r_mat = fi_z_fi_z_r_mat / n_tim;

            a_x_fi_mat = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, init_ang_x_rad_mat, ...
                                            src_frq, cen_frq, lamb, is_eulr_ang);
            a_x_fi_mat = a_x_fi_mat(2: 1: n_sen + 1, :);
            a_z_fi_mat = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, init_ang_z_rad_mat, ...
                                            src_frq, cen_frq, lamb, is_eulr_ang);
            a_z_fi_mat = a_z_fi_mat([1, n_sen + 2: 1: 2 * n_sen], :);

            t_x_mat = [a_x_fc_mat, b_x_fc_mat] * inv([a_x_fi_mat, b_x_fi_mat]);
            t_z_mat = [a_z_fc_mat, b_z_fc_mat] * inv([a_z_fi_mat, b_z_fi_mat]);

            r_x_mat = r_x_mat + (t_x_mat * fi_x_fi_x_r_mat * t_x_mat');
            r_z_mat = r_z_mat + (t_z_mat * fi_z_fi_z_r_mat * t_z_mat');
        end
        
        [eig_vec_x_mat, eig_val_x_mat] = eig(r_x_mat);
        re_eig_val_x_vec               = zeros(n_sen, 1);
        for i_sen = 1: 1: n_sen
            re_eig_val_x_vec(i_sen, 1) = real(eig_val_x_mat(i_sen, i_sen));
        end
        [~, idx_x_vec]     = sort(re_eig_val_x_vec, 'descend');
        sort_eig_vec_x_mat = zeros(n_sen, n_sen);
        for i_sen = 1: 1: n_sen
            sort_eig_vec_x_mat(:, i_sen) = eig_vec_x_mat(:, idx_x_vec(i_sen, 1));
        end
        e_x_mat     = sort_eig_vec_x_mat(:, n_src + 1: n_sen);
        sqr_e_x_mat = e_x_mat * e_x_mat';
        
        [eig_vec_z_mat, eig_val_z_mat] = eig(r_z_mat);
        re_eig_val_z_vec               = zeros(n_sen, 1);
        for i_sen = 1: 1: n_sen
            re_eig_val_z_vec(i_sen, 1) = real(eig_val_z_mat(i_sen, i_sen));
        end
        [~, idx_z_vec]     = sort(re_eig_val_z_vec, 'descend');
        sort_eig_vec_z_mat = zeros(n_sen, n_sen);
        for i_sen = 1: 1: n_sen
            sort_eig_vec_z_mat(:, i_sen) = eig_vec_z_mat(:, idx_z_vec(i_sen, 1));
        end
        e_z_mat     = sort_eig_vec_z_mat(:, n_src + 1: n_sen);
        sqr_e_z_mat = e_z_mat * e_z_mat';
        
        n_ang      = size(ang_rad_vec, 1);
        norm_x_vec = zeros(n_ang, 1);
        norm_z_vec = zeros(n_ang, 1);
        for i_ang = 1: 1: n_ang
            iter_ang_rad_vec = [ang_rad_vec(i_ang, 1), 0, ang_rad_vec(i_ang, 1)];

            ster_vec   = l_fwd_ster_fcn(n_sen_vec, d_sen, iter_ang_rad_vec, lamb, is_eulr_ang);
            ster_x_vec = ster_vec(2: 1: n_sen + 1, :);
            ster_z_vec = ster_vec([1, n_sen + 2: 1: 2 * n_sen], :);

            norm_x_vec(i_ang, 1) = 1 / (ster_x_vec' * sqr_e_x_mat * ster_x_vec);
            norm_z_vec(i_ang, 1) = 1 / (ster_z_vec' * sqr_e_z_mat * ster_z_vec);
        end
        
        [pks_x_vec, locs_x_vec] = findpeaks(abs(norm_x_vec));
        [~,         idx_x_vec]  = sort(pks_x_vec, 'descend');
        [pks_z_vec, locs_z_vec] = findpeaks(abs(norm_z_vec));
        [~,         idx_z_vec]  = sort(pks_z_vec, 'descend');
        
        init_n_src_x = size(pks_x_vec, 1);
        if init_n_src_x > n_src
            init_n_src_x = n_src;
        end
        init_n_src_z = size(pks_z_vec, 1);
        if init_n_src_z > n_src
            init_n_src_z = n_src;
        end
        init_ang_x_rad_vec = ang_rad_vec(locs_x_vec(idx_x_vec(1: init_n_src_x, 1), 1), 1);
        init_ang_z_rad_vec = ang_rad_vec(locs_z_vec(idx_z_vec(1: init_n_src_z, 1), 1), 1);
    end
    
    n_lost_x = n_src - init_n_src_x;
    if n_lost_x > 0
        init_ang_x_rad_vec = [zeros(n_lost_x, 1); init_ang_x_rad_vec];
    end
    n_lost_z = n_src - init_n_src_z;
    if n_lost_z > 0
        init_ang_z_rad_vec = [zeros(n_lost_z, 1); init_ang_z_rad_vec];
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
