%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : ????
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 20 December 2018, Bandhit Suksiri,
%               Updated: 13 March    2019, Bandhit Suksiri.
%
% Copyright 2018 - 2019,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [calc_ang_x_rad_vec, calc_ang_z_rad_vec] = abs_doa_frq_n_src ...
        (src_mat3, frq_vec, pos_mat, lamb, n_src, cen_frq, ang_rad_vec)
    
    [~, idx_ref]          = min(abs(frq_vec - cen_frq));
    [n_frq, n_tim, n_sen] = size(src_mat3);
    n_rank_fc             = n_src;
    n_rank_gsvd           = n_sen;
    src_fc_frq            = frq_vec(idx_ref, 1);
    
    fc_src_mat           = squeeze(src_mat3(idx_ref, :, :)).';
    zero_mean_fc_src_mat = bsxfun(@minus, fc_src_mat, mean(fc_src_mat, 2));
    
    full_w_mat = zeros(n_frq * n_sen,n_sen);
    for i_frq = 1: 1: n_frq
        fi_src_mat           = squeeze(src_mat3(i_frq, :, :)).';
        zero_mean_fi_src_mat = bsxfun(@minus, fi_src_mat, mean(fi_src_mat, 2));
        
        full_w_mat(1 + ((i_frq - 1) * n_sen): i_frq * n_sen, :) = ...
            (zero_mean_fi_src_mat * zero_mean_fc_src_mat') * ...
            (zero_mean_fc_src_mat * zero_mean_fc_src_mat')';
    end
    full_w_mat = full_w_mat ./ n_tim;
    
    v_w_mat   = hgsvd(full_w_mat, n_rank_gsvd, false, true, false);
    sqr_e_mat = v_w_mat(:, end - n_rank_fc: -1: 1) * v_w_mat(:, end - n_rank_fc: -1: 1)';
    
    n_ang    = size(ang_rad_vec, 1);
    norm_mat = zeros(n_ang, n_ang);
    for i_ang_x = 1: 1: n_ang
        for i_ang_z = 1: 1: n_ang
            iter_ang_rad_vec = [ang_rad_vec(i_ang_x, 1), 0, ang_rad_vec(i_ang_z, 1)];
            
            ster_vec = fwd_ster_frq_fcn(pos_mat, iter_ang_rad_vec, src_fc_frq, cen_frq, lamb, true);
            
            norm_mat(i_ang_x, i_ang_z) = (ster_vec' * ster_vec) ./ (ster_vec' * sqr_e_mat * ster_vec);
        end
    end
    abs_norm_mat = abs(norm_mat);
    
    [calc_ang_x_rad_vec, calc_ang_z_rad_vec] = get_peak_2d(abs_norm_mat, ang_rad_vec, n_src);

%     h = surf(abs_norm_mat);
%     set(h, 'edgecolor', 'none');
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%