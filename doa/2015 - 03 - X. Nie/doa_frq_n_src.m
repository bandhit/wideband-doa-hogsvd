%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : X. Nie et al., Array Aperture Extension Algorithm for 2-D DOA Estimation with 
%               L-Shaped Array, 2015.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 10 July 2017, Bandhit Suksiri,
%               Updated: 10 July 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out_ang_x_rad_vec, out_ang_z_rad_vec] = doa_frq_n_src ( ...
    src_mat, n_sen, d_sen, lamb, n_src, src_frq, cen_frq)
    % n_sen_vec must be : [n_sen; 0; n_sen - 1]!
    [n_tim, ~]  = size(src_mat);
    x_mat       = src_mat(:, 2: 1: n_sen + 1).';
    z_mat       = src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
    
    r_mat = zeros(n_sen, n_sen);
    for i = 1: 1: n_tim
        r_mat = r_mat + (x_mat(:, i) * z_mat(:, i)');
    end
    r_mat = r_mat / n_tim;
    
    [u_mat, s_mat, v_mat] = svd(r_mat);
    
    u_1_mat = u_mat(1: n_sen - 1, :);
    u_2_mat = u_mat(2: n_sen,     :);
    v_1_mat = v_mat(1: n_sen - 1, :);
    v_2_mat = v_mat(2: n_sen,     :);
    j_mat   = fliplr(eye(n_sen));
    
    r_1_1_mat = j_mat * ...
                conj(u_mat * pinv(u_1_mat) * u_2_mat * ...
                     s_mat * ...
                     (v_mat * pinv(v_1_mat) * v_2_mat)') * ...
                j_mat;
    r_1_2_mat = u_mat * ((pinv(u_2_mat) * u_1_mat) ^ n_sen) * ...
                s_mat * ...
                v_mat';
    r_2_1_mat = u_mat * ...
                s_mat * ...
                (v_mat * ((pinv(v_2_mat) * v_1_mat) ^ n_sen))';
	
	r_new_mat = [r_1_1_mat, r_1_2_mat; r_2_1_mat, r_mat];

    idx_src_cvec                      = 1: 1: 2 * n_src;
    n_rank                            = size(idx_src_cvec, 2);
    [u_pim_mat, s_pim_mat, v_pim_mat] = svd(r_new_mat);
    
    u_pim_mat = u_pim_mat(:, idx_src_cvec);
    s_pim_mat = s_pim_mat(idx_src_cvec, idx_src_cvec);
    v_pim_mat = v_pim_mat(:, idx_src_cvec);
    new_s_mat = s_pim_mat ^ (1 / 2);
    new_u_mat = u_pim_mat * new_s_mat;
    new_v_mat = v_pim_mat * new_s_mat;
    
    new_u_1_mat = new_u_mat(1: end - 1, :);
    new_u_2_mat = new_u_mat(2: end,     :);
    new_v_1_mat = new_v_mat(1: end - 1, :);
    new_v_2_mat = new_v_mat(2: end,     :);
    
    [p_1_prim_mat, eig_val_u_mat] = eig(pinv(new_u_1_mat) * new_u_2_mat);
    [p_2_prim_mat, eig_val_v_mat] = eig(pinv(new_v_1_mat) * new_v_2_mat);
    eig_val_u_vec                 = diag(eig_val_u_mat);
    eig_val_v_vec                 = diag(eig_val_v_mat);
    
    pair_obj_fcn_1_u_vec = angle(eig_val_u_vec);
    pair_obj_fcn_1_v_vec = angle(eig_val_v_vec);
    pair_u_vec           = zeros(n_src, 2);
    pair_v_vec           = zeros(n_src, 2);
    idx_pair_u           = 0;
    idx_pair_v           = 0;
    for i = 1: 1: n_rank
        [~, min_idx_u_vec] = sort(abs(pair_obj_fcn_1_u_vec(i, 1) - ...
                                      pair_obj_fcn_1_u_vec), 'ascend');
        [~, min_idx_v_vec] = sort(abs(pair_obj_fcn_1_v_vec(i, 1) - ...
                                      pair_obj_fcn_1_v_vec), 'ascend');
        
        is_not_in_pair_u_first_left   = ~any(abs(min_idx_u_vec(1, 1) - pair_u_vec(:, 1)) < eps);
        is_not_in_pair_u_first_right  = ~any(abs(min_idx_u_vec(1, 1) - pair_u_vec(:, 2)) < eps);
        is_not_in_pair_u_second_left  = ~any(abs(min_idx_u_vec(2, 1) - pair_u_vec(:, 1)) < eps);
        is_not_in_pair_u_second_right = ~any(abs(min_idx_u_vec(2, 1) - pair_u_vec(:, 2)) < eps);
        is_not_in_pair_v_first_left   = ~any(abs(min_idx_v_vec(1, 1) - pair_v_vec(:, 1)) < eps);
        is_not_in_pair_v_first_right  = ~any(abs(min_idx_v_vec(1, 1) - pair_v_vec(:, 2)) < eps);
        is_not_in_pair_v_second_left  = ~any(abs(min_idx_v_vec(2, 1) - pair_v_vec(:, 1)) < eps);
        is_not_in_pair_v_second_right = ~any(abs(min_idx_v_vec(2, 1) - pair_v_vec(:, 2)) < eps);
        
        if is_not_in_pair_u_first_left && ...
           is_not_in_pair_u_first_right && ...
           is_not_in_pair_u_second_left && ...
           is_not_in_pair_u_second_right && ...
           (idx_pair_u < n_src)
            idx_pair_u = idx_pair_u + 1;
            pair_u_vec(idx_pair_u, :) = min_idx_u_vec(1: 2, 1).';
        end
        if is_not_in_pair_v_first_left && ...
           is_not_in_pair_v_first_right && ...
           is_not_in_pair_v_second_left && ...
           is_not_in_pair_v_second_right && ...
           (idx_pair_v < n_src)
            idx_pair_v = idx_pair_v + 1;
            pair_v_vec(idx_pair_v, :) = min_idx_v_vec(1: 2, 1).';
        end
        
        if (idx_pair_u == n_src) && (idx_pair_v == n_src)
            break
        end
    end
    p_mat                              = p_1_prim_mat' * p_2_prim_mat;
    raw_pair_obj_fcn_2_mat             = 1 - abs(p_mat);
    tmp_1_pair_obj_fcn_2_mat           = zeros(n_src, n_rank);
    tmp_2_pair_obj_fcn_2_mat           = zeros(n_src, n_src);
    for i = 1: 1: n_src
        tmp_1_pair_obj_fcn_2_mat(i, :) = raw_pair_obj_fcn_2_mat(pair_u_vec(i, 1), :) ...
                                       + raw_pair_obj_fcn_2_mat(pair_u_vec(i, 2), :);
    end
    for i = 1: 1: n_src
        tmp_2_pair_obj_fcn_2_mat(:, i) = tmp_1_pair_obj_fcn_2_mat(:, pair_v_vec(i, 1)) ...
                                       + tmp_1_pair_obj_fcn_2_mat(:, pair_v_vec(i, 1));
    end
    pair_obj_fcn_2_mat                 = tmp_2_pair_obj_fcn_2_mat / 4;
    [~, min_idx_vec]                   = sort(pair_obj_fcn_2_mat(:), 'ascend');
    [row_min_idx_vec, col_min_idx_vec] = ind2sub([n_src, n_src], min_idx_vec);
    sel_row_vec = [];
    sel_col_vec = [];
    for i = 1: 1: n_src * n_src
        row_idx            = row_min_idx_vec(i, 1);
        col_idx            = col_min_idx_vec(i, 1);
        is_not_in_list_row = ~any(abs(row_idx - sel_row_vec) < eps);
        is_not_in_list_col = ~any(abs(col_idx - sel_col_vec) < eps);
        if is_not_in_list_row && ...
           is_not_in_list_col
            sel_row_vec(end + 1, 1) = row_idx; %#ok<AGROW>
            sel_col_vec(end + 1, 1) = col_idx; %#ok<AGROW>
        end
        if (size(sel_row_vec, 1) == n_src) && ...
           (size(sel_col_vec, 1) == n_src)
            break
        end
    end
    
    ang_a_rad_vec = acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                         angle(eig_val_u_vec(pair_u_vec(sel_row_vec, 1), 1)));
	ang_b_rad_vec = acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                         angle(eig_val_v_vec(pair_v_vec(sel_col_vec, 1), 1)));
    ang_a_rad_vec = mean([ang_a_rad_vec, ...
                          acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                               angle(eig_val_u_vec(pair_u_vec(sel_row_vec, 2), 1)))], 2);
	ang_b_rad_vec = mean([ang_b_rad_vec, ...
                          acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                               angle(eig_val_v_vec(pair_v_vec(sel_col_vec, 2), 1)))], 2);
    
    out_ang_x_rad_vec = ang_a_rad_vec;
    out_ang_z_rad_vec = ang_b_rad_vec;
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%