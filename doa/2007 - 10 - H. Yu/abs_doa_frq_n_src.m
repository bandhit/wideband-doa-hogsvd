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
    abs_doa_frq_n_src (src_mat3, frq_vec, pos_mat, lamb, n_src, cen_frq, ang_rad_vec)
    [n_frq, n_tim, n_sen] = size(src_mat3);
    is_eulr_ang           = true;
    
    sqr_e_mat3 = zeros(n_sen, n_sen, n_frq);
    for i_frq = 1: 1: n_frq
        fi_src_mat = squeeze(src_mat3(i_frq, :, :)).';
        fi_x_mat   = fi_src_mat;
        
        fi_x_fi_x_r_mat = zeros(n_sen, n_sen);
        for i_tim = 1: 1: n_tim
            fi_x_fi_x_r_mat = fi_x_fi_x_r_mat + (fi_x_mat(:, i_tim) * fi_x_mat(:, i_tim)');
        end
        fi_x_fi_x_r_mat = fi_x_fi_x_r_mat / n_tim;
        
        r_x_mat = fi_x_fi_x_r_mat;
        
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
        
        sqr_e_mat3(:, :, i_frq) = sqr_e_x_mat;
    end
    
    n_ang    = size(ang_rad_vec, 1);
    norm_mat = zeros(n_ang, n_ang);
    for i_ang_x = 1: 1: n_ang
        for i_ang_z = 1: 1: n_ang
            iter_ang_rad_vec = [ang_rad_vec(i_ang_x, 1), 0, ang_rad_vec(i_ang_z, 1)];
            
            d_vec = zeros(n_frq, 1);
            for i_frq = 1: 1: n_frq
                src_frq   = frq_vec(i_frq, 1);
                sqr_e_mat = sqr_e_mat3(:, :, i_frq);

                ster_vec = fwd_ster_frq_fcn(pos_mat, iter_ang_rad_vec, src_frq, cen_frq, lamb, is_eulr_ang);

                d_vec(i_frq, 1) = ster_vec' * sqr_e_mat * ster_vec;
            end

            norm_mat(i_ang_x, i_ang_z) = 1 / norm(d_vec, 2);
        end
    end
    
    [ang_x_rad_vec, ang_z_rad_vec] = get_peak_2d(abs(norm_mat), ang_rad_vec, n_src);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
