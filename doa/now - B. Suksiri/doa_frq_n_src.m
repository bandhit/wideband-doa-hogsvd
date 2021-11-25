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

function [calc_ang_x_rad_vec, calc_ang_z_rad_vec] = doa_frq_n_src ...
        (src_mat3, frq_vec, n_sen, d_sen, lamb, n_src, cen_frq, doa_mode)
    
    [~, idx_ref]      = min(abs(frq_vec - cen_frq));
    [n_frq, n_tim, ~] = size(src_mat3);
    n_rank_fc         = n_src;
    n_rank_gsvd       = n_sen;
    n_sen_vec         = [n_sen; 0; n_sen - 1];
    src_fc_frq        = frq_vec(idx_ref, 1);
    
    fc_src_mat         = squeeze(src_mat3(idx_ref, :, :)).';
    fc_x_mat           = fc_src_mat(2: 1: n_sen + 1, :);
    fc_z_mat           = fc_src_mat([1, n_sen + 2: 1: 2 * n_sen], :);
    zero_mean_fc_x_mat = bsxfun(@minus, fc_x_mat, mean(fc_x_mat, 2));
    zero_mean_fc_z_mat = bsxfun(@minus, fc_z_mat, mean(fc_z_mat, 2));
    
    full_x_mat = zeros(2 * n_frq * n_sen, n_sen);
    full_z_mat = zeros(2 * n_frq * n_sen, n_sen);
    full_w_mat = zeros(n_frq * (2 * n_sen), 2 * n_sen);
    fi_x_mat3  = zeros(n_sen, n_tim, n_frq);
    fi_z_mat3  = zeros(n_sen, n_tim, n_frq);
    for i_frq = 1: 1: n_frq
        fi_src_mat         = squeeze(src_mat3(i_frq, :, :)).';
        fi_x_mat           = fi_src_mat(2: 1: n_sen + 1, :);
        fi_z_mat           = fi_src_mat([1, n_sen + 2: 1: 2 * n_sen], :);
        zero_mean_fi_x_mat = bsxfun(@minus, fi_x_mat, mean(fi_x_mat, 2));
        zero_mean_fi_z_mat = bsxfun(@minus, fi_z_mat, mean(fi_z_mat, 2));
        
        head_x_idx = 1 + ((i_frq - 1) * n_sen);
        tail_x_idx = i_frq * n_sen;
        head_z_idx = head_x_idx + (n_frq * n_sen);
        tail_z_idx = tail_x_idx + (n_frq * n_sen);
        
        full_x_mat(head_x_idx: tail_x_idx, :) = ...
            (zero_mean_fi_x_mat * zero_mean_fc_z_mat') * ...
            (zero_mean_fc_x_mat * zero_mean_fc_z_mat')';
        full_x_mat(head_z_idx: tail_z_idx, :) = ...
            (zero_mean_fi_z_mat * zero_mean_fc_z_mat') * ...
            (zero_mean_fc_x_mat * zero_mean_fc_z_mat')';
        full_z_mat(head_x_idx: tail_x_idx, :) = ...
            (zero_mean_fi_x_mat * zero_mean_fc_x_mat') * ...
            (zero_mean_fc_z_mat * zero_mean_fc_x_mat')';
        full_z_mat(head_z_idx: tail_z_idx, :) = ...
            (zero_mean_fi_z_mat * zero_mean_fc_x_mat') * ...
            (zero_mean_fc_z_mat * zero_mean_fc_x_mat')';
        
        full_w_mat(1 + ((i_frq - 1) * (2 * n_sen)): i_frq * (2 * n_sen), :) = ...
            (fi_src_mat * fc_src_mat') * ...
            (fc_src_mat * fc_src_mat')';
        
        fi_x_mat3(:, :, i_frq) = fi_x_mat;
        fi_z_mat3(:, :, i_frq) = fi_z_mat;
    end
    full_x_mat = full_x_mat ./ n_tim;
    full_z_mat = full_z_mat ./ n_tim;
    full_w_mat = full_w_mat ./ n_tim;
    
    if doa_mode == "music-no-pair"
        is_qr_schr  = false;
        is_pinv     = true;
        ang_rad_vec = (0: 0.1: 180).' .* pi / 180;
        
        v_x_mat = hgsvd(full_x_mat, n_rank_gsvd, is_qr_schr, is_pinv, false);
        v_z_mat = hgsvd(full_z_mat, n_rank_gsvd, is_qr_schr, is_pinv, false);
        
        norm_x_mat     = music_doa(v_x_mat, n_sen_vec, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, ang_rad_vec, ndef_ang, 2: 1: n_sen + 1);
        abs_norm_x_vec = abs(norm_x_mat(:, 1));
        
        norm_z_mat     = music_doa(v_z_mat, n_sen_vec, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, ndef_ang, ang_rad_vec, [1, n_sen + 2: 1: 2 * n_sen]);
        abs_norm_z_vec = abs(norm_z_mat(1, :)).';
        
        calc_ang_x_rad_vec = get_peak(abs_norm_x_vec, ang_rad_vec, n_rank_fc);
        calc_ang_z_rad_vec = get_peak(abs_norm_z_vec, ang_rad_vec, n_rank_fc);
    elseif doa_mode == "music-signal-auto-pair-norm"
        is_qr_schr  = true;
        is_pinv     = false;
        ang_rad_vec = (0: 1: 180).' .* pi / 180;
        
        v_w_mat = hgsvd(full_w_mat, 2 * n_sen, is_qr_schr, is_pinv, false);
        
        norm_mat     = music_doa_sig(v_w_mat, n_sen_vec, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, ang_rad_vec, ang_rad_vec, 1: 1: 2 * n_sen);
        abs_norm_mat = abs(norm_mat);
        
        [calc_ang_x_rad_vec, calc_ang_z_rad_vec] = get_peak_2d(abs_norm_mat, ang_rad_vec, n_src);
        
        h = surf(abs_norm_mat);
        set(h, 'edgecolor', 'none');
    elseif doa_mode == "music-auto-pair-norm"
        is_qr_schr  = false;
        is_pinv     = true;
        ang_rad_vec = (0: 1: 180).' .* pi / 180;
        
        v_w_mat = hgsvd(full_w_mat, 2 * n_sen, is_qr_schr, is_pinv, false);
        
        norm_mat     = music_doa_sig(v_w_mat, n_sen_vec, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, ang_rad_vec, ang_rad_vec, 1: 1: 2 * n_sen);
        abs_norm_mat = abs(norm_mat);
        
        [calc_ang_x_rad_vec, calc_ang_z_rad_vec] = get_peak_2d(abs_norm_mat, ang_rad_vec, n_src);
        
        h = surf(abs_norm_mat);
        set(h, 'edgecolor', 'none');
    elseif doa_mode == "music-auto-pair-src-var"
        is_qr_schr  = false;
        is_pinv     = true;
        ang_rad_vec = (0: 0.1: 180).' .* pi / 180;
        
        [v_x_mat, u_x_mat3] = hgsvd(full_x_mat, n_rank_gsvd, is_qr_schr, is_pinv, true);
        [v_z_mat, u_z_mat3] = hgsvd(full_z_mat, n_rank_gsvd, is_qr_schr, is_pinv, true);
        
        norm_x_mat     = music_doa(v_x_mat, n_sen_vec, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, ang_rad_vec, ndef_ang, 2: 1: n_sen + 1);
        abs_norm_x_vec = abs(norm_x_mat(:, 1));
        
        norm_z_mat     = music_doa(v_z_mat, n_sen_vec, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, ndef_ang, ang_rad_vec, [1, n_sen + 2: 1: 2 * n_sen]);
        abs_norm_z_vec = abs(norm_z_mat(1, :)).';
        
        calc_ang_x_rad_vec = get_peak(abs_norm_x_vec, ang_rad_vec, n_rank_fc);
        calc_ang_z_rad_vec = get_peak(abs_norm_z_vec, ang_rad_vec, n_rank_fc);
        
        ang_rad_mat = [calc_ang_x_rad_vec, ndef_ang .* ones(n_rank_fc, 1), calc_ang_z_rad_vec];
        
        a_fi_mat   = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, ang_rad_mat, src_fc_frq, cen_frq, lamb, true);
        a_fi_x_mat = a_fi_mat(2: 1: n_sen + 1, :);
        a_fi_z_mat = a_fi_mat([1, n_sen + 2: 1: 2 * n_sen], :);
        
        u_x_fc_x_mat          = u_x_mat3(:, end - n_rank_fc + 1: end, idx_ref);
        u_z_fc_z_mat          = u_z_mat3(:, end - n_rank_fc + 1: end, idx_ref + n_frq);
        src_fi_z_src_fi_z_mat = zeros(n_src, n_src);
        src_fi_x_src_fi_z_mat = zeros(n_src, n_src);
        for i_frq = 1: 1: n_frq
            u_x_fi_x_mat = u_x_mat3(:, end - n_rank_fc + 1: end, i_frq);
            u_z_fi_z_mat = u_z_mat3(:, end - n_rank_fc + 1: end, i_frq + n_frq);
            
            t_x_mat = u_x_fc_x_mat / u_x_fi_x_mat;
            t_z_mat = u_z_fc_z_mat / u_z_fi_z_mat;
            
            fi_x_mat = fi_x_mat3(:, :, i_frq);
            fi_z_mat = fi_z_mat3(:, :, i_frq);
            
            src_fi_x_mat = a_fi_x_mat \ t_x_mat * fi_x_mat;
            src_fi_z_mat = a_fi_z_mat \ t_z_mat * fi_z_mat;
            
            src_fi_z_src_fi_z_mat = src_fi_z_src_fi_z_mat + (src_fi_z_mat * src_fi_z_mat');
            src_fi_x_src_fi_z_mat = src_fi_x_src_fi_z_mat + (src_fi_x_mat * src_fi_z_mat');
        end
        src_fi_z_src_fi_z_mat = src_fi_z_src_fi_z_mat ./ n_tim;
        src_fi_x_src_fi_z_mat = src_fi_x_src_fi_z_mat ./ n_tim;
        opp_mat               = src_fi_z_src_fi_z_mat * src_fi_x_src_fi_z_mat';
        
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

        calc_ang_x_rad_vec = calc_ang_x_rad_vec(sel_col_fc_vec(sel_row_fc_vec), 1);
    elseif doa_mode == "esprit-no-pair"
        is_qr_schr = true;
        is_pinv    = false;
        
        v_x_mat = hgsvd(full_x_mat, n_rank_gsvd, is_qr_schr, is_pinv, false);
        v_z_mat = hgsvd(full_z_mat, n_rank_gsvd, is_qr_schr, is_pinv, false);
        
        calc_ang_x_rad_vec = esprit_doa(v_x_mat, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, true, true);
        calc_ang_z_rad_vec = esprit_doa(v_z_mat, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, true, true);
    elseif doa_mode == "esprit-auto-pair-x"
        is_qr_schr = true;
        is_pinv    = false;
        
        [v_x_mat, u_x_mat3] = hgsvd(full_x_mat, n_rank_gsvd, is_qr_schr, is_pinv, true);
        u_x_fc_x_mat        = u_x_mat3(:, end - n_rank_fc + 1: end, idx_ref);
        u_x_fc_z_mat        = u_x_mat3(:, end - n_rank_fc + 1: end, idx_ref + n_frq);
        
        calc_ang_x_rad_vec = esprit_doa(v_x_mat, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, true, true);
        
        ang_rad_mat = [calc_ang_x_rad_vec, ndef_ang .* ones(n_rank_fc, 1), ndef_ang .* ones(n_rank_fc, 1)];
        
        a_fi_mat   = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, ang_rad_mat, src_fc_frq, cen_frq, lamb, true);
        a_fi_x_mat = a_fi_mat(2: 1: n_sen + 1, :);
        
        a_fi_z_mat = u_x_fc_z_mat / u_x_fc_x_mat * a_fi_x_mat;
        phi_z_vec  = diag(lsqminnorm(a_fi_z_mat(1: end - 1, :), a_fi_z_mat(2: end, :)));
        
        calc_ang_z_rad_vec = acos(((cen_frq * lamb) / (2 * pi * d_sen * src_fc_frq)) * angle(phi_z_vec));
    elseif doa_mode == "esprit-auto-pair-z"
        is_qr_schr = true;
        is_pinv    = false;
        
        [v_z_mat, u_z_mat3] = hgsvd(full_z_mat, n_rank_gsvd, is_qr_schr, is_pinv, true);
        u_z_fc_z_mat        = u_z_mat3(:, end - n_rank_fc + 1: end, idx_ref + n_frq);
        u_z_fc_x_mat        = u_z_mat3(:, end - n_rank_fc + 1: end, idx_ref);
        
        calc_ang_z_rad_vec = esprit_doa(v_z_mat, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, true, true);
        
        ang_rad_mat = [ndef_ang .* ones(n_rank_fc, 1), ndef_ang .* ones(n_rank_fc, 1), calc_ang_z_rad_vec];
        
        a_fi_mat   = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, ang_rad_mat, src_fc_frq, cen_frq, lamb, true);
        a_fi_z_mat = a_fi_mat([1, n_sen + 2: 1: 2 * n_sen], :);
        
        a_fi_x_mat = u_z_fc_x_mat / u_z_fc_z_mat * a_fi_z_mat;
        phi_x_vec  = diag(lsqminnorm(a_fi_x_mat(1: end - 1, :), a_fi_x_mat(2: end, :)));
        
        calc_ang_x_rad_vec = acos(((cen_frq * lamb) / (2 * pi * d_sen * src_fc_frq)) * angle(phi_x_vec));
    elseif doa_mode == "esprit-auto-pair-src-var"
        is_qr_schr = true;
        is_pinv    = false;
        
        [v_x_mat, u_x_mat3] = hgsvd(full_x_mat, n_rank_gsvd, is_qr_schr, is_pinv, true);
        [v_z_mat, u_z_mat3] = hgsvd(full_z_mat, n_rank_gsvd, is_qr_schr, is_pinv, true);
        
%         trace(u_x_mat3(:, :, 100) * u_x_mat3(:, :, 100)')
%         trace(u_x_mat3(:, :, 100)' * u_x_mat3(:, :, 100))
%         trace(u_x_mat3(:, :, 100)  \ u_x_mat3(:, :, 100))
%         trace(inv(u_x_mat3(:, :, 100) * u_x_mat3(:, :, 100)'))
%         trace(inv(u_x_mat3(:, :, 100)' * u_x_mat3(:, :, 100)))
        
        calc_ang_x_rad_vec = esprit_doa(v_x_mat, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, true, true);
        calc_ang_z_rad_vec = esprit_doa(v_z_mat, src_fc_frq, cen_frq, lamb, d_sen, n_rank_fc, true, true);
        
        ang_rad_mat = [calc_ang_x_rad_vec, ndef_ang .* ones(n_rank_fc, 1), calc_ang_z_rad_vec];
        
        a_fi_mat   = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, ang_rad_mat, src_fc_frq, cen_frq, lamb, true);
        a_fi_x_mat = a_fi_mat(2: 1: n_sen + 1, :);
        a_fi_z_mat = a_fi_mat([1, n_sen + 2: 1: 2 * n_sen], :);
        
%         (a_fi_x_mat' * a_fi_x_mat)
%         trace((a_fi_x_mat' * a_fi_x_mat) * (a_fi_x_mat' * a_fi_x_mat))
        
        u_x_fc_x_mat          = u_x_mat3(:, end - n_rank_fc + 1: end, idx_ref);
        u_z_fc_z_mat          = u_z_mat3(:, end - n_rank_fc + 1: end, idx_ref + n_frq);
        src_fi_z_src_fi_z_mat = zeros(n_src, n_src);
        src_fi_x_src_fi_z_mat = zeros(n_src, n_src);
        for i_frq = 1: 1: n_frq
            u_x_fi_x_mat = u_x_mat3(:, end - n_rank_fc + 1: end, i_frq);
            u_z_fi_z_mat = u_z_mat3(:, end - n_rank_fc + 1: end, i_frq + n_frq);
            
            t_x_mat = u_x_fc_x_mat / u_x_fi_x_mat;
            t_z_mat = u_z_fc_z_mat / u_z_fi_z_mat;
            
            fi_x_mat = fi_x_mat3(:, :, i_frq);
            fi_z_mat = fi_z_mat3(:, :, i_frq);
            
            src_fi_x_mat = a_fi_x_mat \ t_x_mat * fi_x_mat;
            src_fi_z_mat = a_fi_z_mat \ t_z_mat * fi_z_mat;
            
            src_fi_z_src_fi_z_mat = src_fi_z_src_fi_z_mat + (src_fi_z_mat * src_fi_z_mat');
            src_fi_x_src_fi_z_mat = src_fi_x_src_fi_z_mat + (src_fi_x_mat * src_fi_z_mat');
        end
        src_fi_z_src_fi_z_mat = src_fi_z_src_fi_z_mat ./ n_tim;
        src_fi_x_src_fi_z_mat = src_fi_x_src_fi_z_mat ./ n_tim;
        opp_mat               = src_fi_z_src_fi_z_mat * src_fi_x_src_fi_z_mat';
        
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

        calc_ang_x_rad_vec = calc_ang_x_rad_vec(sel_col_fc_vec(sel_row_fc_vec), 1);
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%