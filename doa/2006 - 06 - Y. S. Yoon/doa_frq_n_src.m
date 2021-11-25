%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : Y. S. Yoon et al., TOPS: New DOA Estimator for Wideband Signals, IEEE, 2006.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 11 September 2017, Bandhit Suksiri,
%               Updated: 12 September 2017, Bandhit Suksiri.
%
% Copyright 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ang_x_rad_vec, ang_z_rad_vec] = ...
    doa_frq_n_src (src_mat3, frq_vec, n_sen, d_sen, lamb, n_src, cen_frq, idx_ref, ang_rad_vec)
    [n_frq, n_tim, ~] = size(src_mat3);
    
    f_x_mat3 = zeros(n_sen, n_src,         n_frq);
    f_z_mat3 = zeros(n_sen, n_src,         n_frq);
    w_x_mat3 = zeros(n_sen, n_sen - n_src, n_frq);
    w_z_mat3 = zeros(n_sen, n_sen - n_src, n_frq);
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
        f_x_mat = sort_eig_vec_x_mat(:,         1: n_src);
        w_x_mat = sort_eig_vec_x_mat(:, n_src + 1: n_sen);
        
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
        f_z_mat = sort_eig_vec_z_mat(:,         1: n_src);
        w_z_mat = sort_eig_vec_z_mat(:, n_src + 1: n_sen);
        
        f_x_mat3(:, :, i_frq) = f_x_mat;
        w_x_mat3(:, :, i_frq) = w_x_mat;
        f_z_mat3(:, :, i_frq) = f_z_mat;
        w_z_mat3(:, :, i_frq) = w_z_mat;
    end
    
    src_frq_ref = frq_vec(idx_ref, 1);
    f_x_ref_mat = f_x_mat3(:, :, idx_ref);
    f_z_ref_mat = f_z_mat3(:, :, idx_ref);
    
    m_x_vec = (1: 1: n_sen).';
    m_z_vec = (0: 1: n_sen - 1).';
    
    n_ang      = size(ang_rad_vec, 1);
    norm_x_vec = zeros(n_ang, 1);
    norm_z_vec = zeros(n_ang, 1);
    for i_ang = 1: 1: n_ang
        d_x_mat = [];
        d_z_mat = [];
        for i_frq = 1: 1: n_frq
            src_frq  = frq_vec(i_frq, 1);
            diff_frq = abs(src_frq - src_frq_ref);
            w_x_mat  = w_x_mat3(:, :, i_frq);
            w_z_mat  = w_z_mat3(:, :, i_frq);
            
            ohm_x_mat = frq_tran_mat( ...
                m_x_vec, d_sen, diff_frq, cen_frq, lamb, ang_rad_vec(i_ang, 1));
            ohm_z_mat = frq_tran_mat( ...
                m_z_vec, d_sen, diff_frq, cen_frq, lamb, ang_rad_vec(i_ang, 1));
            
            u_x_mat = ohm_x_mat * f_x_ref_mat;
            z_x_mat = ohm_z_mat * f_z_ref_mat;
            d_x_mat = [d_x_mat, u_x_mat' * w_x_mat]; %#ok<AGROW>
            d_z_mat = [d_z_mat, z_x_mat' * w_z_mat]; %#ok<AGROW>
        end
        norm_x_vec(i_ang, 1) = 1 / min(svd(d_x_mat));
        norm_z_vec(i_ang, 1) = 1 / min(svd(d_z_mat));
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

function ohm_mat = frq_tran_mat (m_vec, d_sen, src_frq, cen_frq, lamb, ang_rad)
    ohm_vec = exp((- (1i * 2 * pi * d_sen * src_frq * cos(ang_rad)) / (cen_frq * lamb)) .* m_vec);
    ohm_mat = diag(ohm_vec);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
