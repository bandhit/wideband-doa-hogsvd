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
% Logs        : Created: 13 June    2017, Bandhit Suksiri,
%               Updated: 11 January 2019, Bandhit Suksiri.
%
% Copyright 2016 - 2019,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ang_x_rad_vec, ang_z_rad_vec] = doa_frq_n_src (src_mat3, frq_vec, n_sen, d_sen, lamb, n_src, cen_frq)
    [~, idx_ref]      = min(abs(frq_vec - cen_frq));
    [n_frq, n_tim, ~] = size(src_mat3);
    n_rank_fi         = n_src;
    n_rank_fc         = n_src;
    src_fc_frq        = frq_vec(idx_ref, 1);
    n_sen_vec         = [n_sen; 0; n_sen - 1];
    
    fc_src_mat = squeeze(src_mat3(idx_ref, :, :));
    fc_x_mat   = fc_src_mat(:, 2: 1: n_sen + 1).';
    fc_z_mat   = fc_src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
    fc_x_mat   = bsxfun(@minus, fc_x_mat, mean(fc_x_mat, 2));
    fc_z_mat   = bsxfun(@minus, fc_z_mat, mean(fc_z_mat, 2));
    
    
    
    
    sum_tran_fi_x_mat = zeros(n_sen, n_tim);
    sum_tran_fi_z_mat = zeros(n_sen, n_tim);
    for i_frq = 1: 1: n_frq
        fi_src_mat = squeeze(src_mat3(i_frq, :, :));
        fi_x_mat   = fi_src_mat(:, 2: 1: n_sen + 1).';
        fi_z_mat   = fi_src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
        fi_x_mat   = bsxfun(@minus, fi_x_mat, mean(fi_x_mat, 2));
        fi_z_mat   = bsxfun(@minus, fi_z_mat, mean(fi_z_mat, 2));
        
        r_x_mat = zeros(n_sen, n_sen);
        r_z_mat = zeros(n_sen, n_sen);
        for i_tim = 1: 1: n_tim
            r_x_mat = r_x_mat ...
                    + (fi_x_mat(:, i_tim)  * ...
                       fc_z_mat(:, i_tim)' * ...
                       fc_z_mat(:, i_tim)  * ...
                       fc_x_mat(:, i_tim)');
            r_z_mat = r_z_mat ...
                    + (fi_z_mat(:, i_tim)  * ...
                       fc_x_mat(:, i_tim)' * ...
                       fc_x_mat(:, i_tim)  * ...
                       fc_z_mat(:, i_tim)');
        end
        
        [r_x_u_bar_fi_mat, ~, r_x_v_bar_fc_mat] = svd(r_x_mat);
        [r_z_u_bar_fi_mat, ~, r_z_v_bar_fc_mat] = svd(r_z_mat);
        
        tran_fi_x_mat = r_x_v_bar_fc_mat(:, 1: n_rank_fi) * r_x_u_bar_fi_mat(:, 1: n_rank_fi)' * fi_x_mat;
        tran_fi_z_mat = r_z_v_bar_fc_mat(:, 1: n_rank_fi) * r_z_u_bar_fi_mat(:, 1: n_rank_fi)' * fi_z_mat;
        
%         tran_fi_x_mat = bsxfun(@minus, tran_fi_x_mat, mean(tran_fi_x_mat, 2));
%         tran_fi_z_mat = bsxfun(@minus, tran_fi_z_mat, mean(tran_fi_z_mat, 2));
        
        sum_tran_fi_x_mat = sum_tran_fi_x_mat + tran_fi_x_mat;
        sum_tran_fi_z_mat = sum_tran_fi_z_mat + tran_fi_z_mat;
    end
    sum_tran_fi_x_mat = bsxfun(@minus, sum_tran_fi_x_mat, mean(sum_tran_fi_x_mat, 2));
    sum_tran_fi_z_mat = bsxfun(@minus, sum_tran_fi_z_mat, mean(sum_tran_fi_z_mat, 2));
    
    
    
    
    [u_r_mat, ~, v_r_mat] = svd(sum_tran_fi_x_mat * sum_tran_fi_z_mat');
    
    u_1_mat = u_r_mat(1: end - 1, 1: n_rank_fc);
    u_2_mat = u_r_mat(2: end,     1: n_rank_fc);
    v_1_mat = v_r_mat(1: end - 1, 1: n_rank_fc);
    v_2_mat = v_r_mat(2: end,     1: n_rank_fc);
    
    eig_val_u_vec = eig(pinv(u_1_mat) * u_2_mat);
    eig_val_v_vec = eig(pinv(v_1_mat) * v_2_mat);
    
    ang_x_rad_vec = acos(((cen_frq * lamb) / (2 * pi * d_sen * src_fc_frq)) * angle(eig_val_u_vec));
	ang_z_rad_vec = acos(((cen_frq * lamb) / (2 * pi * d_sen * src_fc_frq)) * angle(eig_val_v_vec));
    ang_rad_mat   = [ang_x_rad_vec, ndef_ang .* ones(n_rank_fc, 1), ang_z_rad_vec];
    
    
    
    
    a_fi_mat   = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, ang_rad_mat, src_fc_frq, cen_frq, lamb, true);
    a_fi_x_mat = a_fi_mat(2: 1: n_sen + 1, :);
    a_fi_z_mat = a_fi_mat([1, n_sen + 2: 1: 2 * n_sen], :);
    
    src_sum_tran_fi_x_mat = pinv(a_fi_x_mat) * sum_tran_fi_x_mat;
    src_sum_tran_fi_z_mat = pinv(a_fi_z_mat) * sum_tran_fi_z_mat;
    
    src_fi_z_src_fi_z_mat = src_sum_tran_fi_z_mat * src_sum_tran_fi_z_mat';
    src_fi_x_src_fi_z_mat = src_sum_tran_fi_x_mat * src_sum_tran_fi_z_mat';
    
    opp_mat                     = src_fi_z_src_fi_z_mat * src_fi_x_src_fi_z_mat';
    [l_pair_mat, ~, r_pair_mat] = svd(opp_mat);
    pair_mat                    = 1 - abs(l_pair_mat * r_pair_mat');
    
    
    
    
    [~, min_fc_idx_vec]                      = sort(pair_mat(:), 'ascend');
    [row_min_fc_idx_vec, col_min_fc_idx_vec] = ind2sub([n_rank_fc, n_rank_fc], min_fc_idx_vec);
    
    n_fnd_src      = 0;
    sel_row_fc_vec = [];
    sel_col_fc_vec = [];
    for i = 1: 1: n_rank_fc * n_rank_fc
        row_idx            = row_min_fc_idx_vec(i, 1);
        col_idx            = col_min_fc_idx_vec(i, 1);
        is_not_in_list_row = ~any(abs(row_idx - sel_row_fc_vec) < eps);
        is_not_in_list_col = ~any(abs(col_idx - sel_col_fc_vec) < eps);
        if is_not_in_list_row && is_not_in_list_col
            sel_row_fc_vec(end + 1, 1) = row_idx; %#ok<AGROW>
            sel_col_fc_vec(end + 1, 1) = col_idx; %#ok<AGROW>
            n_fnd_src                  = n_fnd_src + 1;
            if n_fnd_src == n_rank_fc
                break
            end
        end
    end
    
    ang_x_rad_vec = ang_x_rad_vec(sel_col_fc_vec(sel_row_fc_vec), 1);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
