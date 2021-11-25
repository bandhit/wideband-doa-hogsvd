%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : H. Yu et al., A New Method for Wideband DOA Estimation, IEEE, 2007.
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

function [ang_x_rad_vec, ang_z_rad_vec] = ...
    doa_frq_n_src (src_mat3, frq_vec, n_sen, d_sen, lamb, n_src, cen_frq, ang_rad_vec)
    [n_frq, n_tim, ~] = size(src_mat3);
    n_sen_vec         = [n_sen; 0; n_sen - 1];
    is_eulr_ang       = true;
    
    sqr_e_x_mat3 = zeros(n_sen, n_sen, n_frq);
    sqr_e_z_mat3 = zeros(n_sen, n_sen, n_frq);
    for i_frq = 1: 1: n_frq
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
        
        r_x_mat = fi_x_fi_x_r_mat;
        r_z_mat = fi_z_fi_z_r_mat;
        
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
        
        sqr_e_x_mat3(:, :, i_frq) = sqr_e_x_mat;
        sqr_e_z_mat3(:, :, i_frq) = sqr_e_z_mat;
    end
    
    n_ang      = size(ang_rad_vec, 1);
    norm_x_vec = zeros(n_ang, 1);
    norm_z_vec = zeros(n_ang, 1);
    for i_ang = 1: 1: n_ang
        iter_ang_rad_vec = [ang_rad_vec(i_ang, 1), 0, ang_rad_vec(i_ang, 1)];
        
        d_x_vec = zeros(n_frq, 1);
        d_z_vec = zeros(n_frq, 1);
        for i_frq = 1: 1: n_frq
            src_frq     = frq_vec(i_frq, 1);
            sqr_e_x_mat = sqr_e_x_mat3(:, :, i_frq);
            sqr_e_z_mat = sqr_e_z_mat3(:, :, i_frq);
            
            ster_vec    = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, iter_ang_rad_vec, ...
                                             src_frq, cen_frq, lamb, is_eulr_ang);
            ster_x_vec  = ster_vec(2: 1: n_sen + 1, :);
            ster_z_vec  = ster_vec([1, n_sen + 2: 1: 2 * n_sen], :);

            d_x_vec(i_frq, 1) = ster_x_vec' * sqr_e_x_mat * ster_x_vec;
            d_z_vec(i_frq, 1) = ster_z_vec' * sqr_e_z_mat * ster_z_vec;
        end
        
        norm_x_vec(i_ang, 1) = 1 / norm(d_x_vec, 2);
        norm_z_vec(i_ang, 1) = 1 / norm(d_z_vec, 2);
    end
    
    [pks_x_vec, locs_x_vec] = findpeaks(abs(norm_x_vec));
    [~,         idx_x_vec]  = sort(pks_x_vec, 'descend');
    [pks_z_vec, locs_z_vec] = findpeaks(abs(norm_z_vec));
    [~,         idx_z_vec]  = sort(pks_z_vec, 'descend');

    post_n_src_x = size(pks_x_vec, 1);
    if post_n_src_x > n_src
        post_n_src_x = n_src;
    end
    post_n_src_z = size(pks_z_vec, 1);
    if post_n_src_z > n_src
        post_n_src_z = n_src;
    end
    ang_x_rad_vec = ang_rad_vec(locs_x_vec(idx_x_vec(1: post_n_src_x, 1), 1), 1);
    ang_z_rad_vec = ang_rad_vec(locs_z_vec(idx_z_vec(1: post_n_src_z, 1), 1), 1);
    
    n_lost_x = n_src - post_n_src_x;
    if n_lost_x > 0
        ang_x_rad_vec = [zeros(n_lost_x, 1); ang_x_rad_vec];
    end
    n_lost_z = n_src - post_n_src_z;
    if n_lost_z > 0
        ang_z_rad_vec = [zeros(n_lost_z, 1); ang_z_rad_vec];
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
