%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : A Computationally Efficient Wideband Direction-of-Arrival Estimation Method for 
%               L-Shaped Microphone Arrays, IEEE, 2018.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 13 June 2017, Bandhit Suksiri,
%               Updated: 16 May  2018, Bandhit Suksiri.
%
% Copyright 2016 - 2018,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ang_x_rad_vec, ang_z_rad_vec, s_r_mat, u_r_mat, v_r_mat] = ...
    doa_frq_n_src (src_mat3, frq_vec, n_sen, d_sen, lamb, n_src, cen_frq, idx_ref)
    [n_frq, n_tim, ~] = size(src_mat3);
    n_rank_fc         = n_src;
    src_frq           = frq_vec(idx_ref, 1);
    fc_src_mat        = squeeze(src_mat3(idx_ref, :, :));
    
    % n_sen_vec must be : [n_sen; 0; n_sen - 1]!
    fc_x_mat = fc_src_mat(:, 2: 1: n_sen + 1).';
    fc_z_mat = fc_src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
    fc_z_fc_z_r_mat = zeros(n_sen, n_sen);
    for i_tim = 1: 1: n_tim
        fc_z_fc_z_r_mat = fc_z_fc_z_r_mat + (fc_z_mat(:, i_tim) * fc_z_mat(:, i_tim)');
    end
    fc_z_fc_z_r_mat    = fc_z_fc_z_r_mat / n_tim;    
    re_fc_z_fc_z_r_mat = fc_z_fc_z_r_mat + fc_z_fc_z_r_mat';
    
    [eig_vec_tmp_z_mat, eig_val_tmp_z_mat] = eig(re_fc_z_fc_z_r_mat);
    [~, idx_tmp_z_vec]                     = sort(real(diag(eig_val_tmp_z_mat)), 'descend');
    s_eig_vec_tmp_z_mat                    = eig_vec_tmp_z_mat(:, idx_tmp_z_vec);
    tmp_r_z_mat                            = s_eig_vec_tmp_z_mat(:, 1: n_rank_fc) ...
                                           * s_eig_vec_tmp_z_mat(:, 1: n_rank_fc)';
    
    r_mat = zeros(n_sen, n_sen);
    for i_frq = 1: 1: n_frq
        fi_src_mat = squeeze(src_mat3(i_frq, :, :));
        fi_x_mat   = fi_src_mat(:, 2: 1: n_sen + 1).';
        fi_z_mat   = fi_src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
        
        % conventional way
        fc_x_fi_z_r_mat = zeros(n_sen, n_sen);
        fi_z_fi_x_r_mat = zeros(n_sen, n_sen);
        fc_z_fi_x_r_mat = zeros(n_sen, n_sen);
        for i_tim = 1: 1: n_tim
            fc_x_fi_z_r_mat = fc_x_fi_z_r_mat + (fc_x_mat(:, i_tim) * fi_z_mat(:, i_tim)');
            fi_z_fi_x_r_mat = fi_z_fi_x_r_mat + (fi_z_mat(:, i_tim) * fi_x_mat(:, i_tim)');
            fc_z_fi_x_r_mat = fc_z_fi_x_r_mat + (fc_z_mat(:, i_tim) * fi_x_mat(:, i_tim)');
            
        end
        fc_x_fi_z_r_mat = fc_x_fi_z_r_mat / n_tim;
        fi_z_fi_x_r_mat = fi_z_fi_x_r_mat / n_tim;
        fc_z_fi_x_r_mat = fc_z_fi_x_r_mat / n_tim;
        tmp_r_mat       = fc_x_fi_z_r_mat * fi_z_fi_x_r_mat / fc_z_fi_x_r_mat;
        
        r_mat = r_mat + tmp_r_mat;
    end
    r_mat = r_mat * tmp_r_z_mat;
    
    [u_r_mat, s_r_mat, v_r_mat] = svd(r_mat);
    s_r_vec                     = diag(s_r_mat);
    
    sqrt_s_r_mat = diag(sqrt(s_r_vec(1: n_rank_fc, 1)));
    
    u_1_mat = u_r_mat(1: end - 1, 1: n_rank_fc) * sqrt_s_r_mat;
    u_2_mat = u_r_mat(2: end,     1: n_rank_fc) * sqrt_s_r_mat;
    v_1_mat = v_r_mat(1: end - 1, 1: n_rank_fc) * sqrt_s_r_mat;
    v_2_mat = v_r_mat(2: end,     1: n_rank_fc) * sqrt_s_r_mat;
    
    [p_1_prim_mat, eig_val_u_mat] = eig(pinv(u_1_mat) * u_2_mat);
    [p_2_prim_mat, eig_val_v_mat] = eig(pinv(v_1_mat) * v_2_mat);
    
    eig_val_u_vec = diag(eig_val_u_mat);
    eig_val_v_vec = diag(eig_val_v_mat);
    
    pair_fc_mat                              = 1 - abs(p_2_prim_mat' * p_1_prim_mat);
    [~, min_fc_idx_vec]                      = sort(pair_fc_mat(:), 'ascend');
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
    
    ang_x_rad_vec = acos(((cen_frq * lamb) / (2 * pi * d_sen * src_frq)) * ...
                         real(log(eig_val_u_vec(sel_col_fc_vec, 1)) ./ 1i));
	ang_z_rad_vec = acos(((cen_frq * lamb) / (2 * pi * d_sen * src_frq)) * ...
                         real(log(eig_val_v_vec(sel_col_fc_vec, 1)) ./ 1i));
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
